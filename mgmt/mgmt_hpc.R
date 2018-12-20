# This script runs simulations of different management plans for the University
# of New Hampshire woodlands, which include 22 properties totaling nearly 3,500
# acres. The code is written to run iterations in parallel on a high performance
# cluster, and the number of cores to use can be specified with
# set_g_p(n.cores=...). We explore four different broad plans:
# 1) No management
# 2) Current stated management plan (property-specific)
# 3) Current practiced management plan (property-specific)
# 3) Manual removal for any populations > 500 individuals/20-acres

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=T,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")

########
## Setup
########
# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW", "fastmatch")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

# set parameters
res <- c("20ac", "9km2")[2]
n.sim <- 12
dem_par <- set_g_p(tmax=20, n.cores=4)
if(res == "9km2") {
  dem_par$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  dem_par$sdd.max <- 7
  dem_par$sdd.rate <- 1.4
}

# load landscape
lc.df <- read.csv(paste0("data/USDA_", res, "_mgmt.csv"))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initial populations
N_0 <- readRDS(paste0("data/inits/N_2018_", res, ".rds"))[lc.df$id.full,,]
B_0 <- readRDS(paste0("data/inits/B_2018_", res, ".rds"))[lc.df$id.full]
sdd <- sdd_set_probs(ncell, lc.df, dem_par, verbose=T)

# management plans
mgmt_plans <- dir("mgmt", "^p", full.names=T) %>% map(read.csv) %>%
  setNames(str_remove(str_remove(dir("mgmt", "^p"), "p._"), ".csv"))
mgmt <- setNames(vector("list", length(mgmt_plans)), names(mgmt_plans))
lc.UNH <- filter(lc.df, in.UNH)
for(m in seq_along(mgmt)) {
  if(names(mgmt)[m] == "none") {
    mgmt[[m]] <- list(ctrl_par=set_control_p(null_ctrl=TRUE),
                      manual.i=NULL,
                      cover.i=NULL)
  } else {
    mgmt_index <- match(lc.UNH$Property, mgmt_plans[[m]]$Property)
    mgmt_m <- cbind(id=lc.UNH$id, mgmt_plans[[m]][mgmt_index,])
    mgmt[[m]] <- list(ctrl_par=set_control_p(null_ctrl=FALSE),
                      manual.i=mgmt_m[,grep("id|Trt", names(mgmt_m))],
                      thresh=mgmt_m$thresh,
                      trt.int=mgmt_m$trt.cycle,
                      trt.followup=mgmt_m$trt.followup,
                      cover.i=mgmt_m[,grep("id|ground", names(mgmt_m))] %>%
                        dplyr::rename(Trt=ground_cover))
  }
} 

# storage and initialization objects
out.df.ls <- out.all.ls <- setNames(vector("list", length(mgmt)), names(mgmt))
N <- array(0, dim=c(ngrid, dem_par$tmax+1, 6, max(dem_par$m)))
N[,1,,] <- N_0
B <- matrix(0, nrow=ngrid, ncol=dem_par$tmax+1)
B[,1] <- B_0

for(m in seq_along(mgmt)) {
  # ids of all managed cells
  managed <- unique(c(mgmt[[m]]$manual.i$id, mgmt[[m]]$cover.i$id))
  
  # storage for simulation output
  N_m <- B_m <- vector("list", n.sim)
  
  # iterate through simulations
  p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
  out_m <- foreach(s=1:n.sim, 
                   .packages=c("gbPopMod", "tidyverse", "magrittr")) %dopar% {
    # iterate through years
    N_s <- N
    B_s <- B
    for(k in 1:dem_par$tmax) {
      # decide which cells are managed this year
      if(names(mgmt)[m] != "none" && 
         any((k %% mgmt[[m]]$trt.int) <= mgmt[[m]]$trt.followup)) {  
        managed_k <- managed[apply(N_s[managed,k,,], 1, sum) > mgmt[[m]]$thresh &
                               (k %% mgmt[[m]]$trt.int) <= mgmt[[m]]$trt.followup]
        manual_k <- mgmt[[m]]$manual.i[mgmt[[m]]$manual.i$id %in% managed_k,]
        cover_k <- mgmt[[m]]$cover.i[mgmt[[m]]$cover.i$id %in% managed_k,]
        if(names(mgmt)[m] != "aggressive") cover_k <- NULL
      } else {
        manual_k <- NULL
        cover_k <- NULL
      }
      
      # buckthorn management, population growth, dispersal
      out <- iterate_pop(ngrid, ncell, N_s[,k,,], B_s[,k], dem_par, 
                         lc.df[,-(15:18)], sdd, mgmt[[m]]$ctrl_par, 
                         cover_k, manual_k)
      
      # update abundances
      N_s[,k+1,,] <- out$N
      B_s[,k+1] <- out$B
    }
    # return abundances
    list(N=N_s, B=B_s)
  }
  stopCluster(p.c)
  N_m <- map(out_m, ~.$N)
  B_m <- map(out_m, ~.$B)
  
  # summarise: means across simulations
  N_adult.tot <- Reduce('+', map(N_m, ~apply(.[,,,max(dem_par$m)], 1:2, sum)))/n.sim
  N_all.tot <- Reduce('+', map(N_m, ~apply(., 1:2, sum)))/n.sim
  B.tot <- Reduce('+', B_m)/n.sim
  N_adult.df <- as.data.frame(N_adult.tot); names(N_adult.df) <- 1:ncol(N_adult.df)
  N_all.df <- as.data.frame(N_all.tot); names(N_all.df) <- 1:ncol(N_all.df)
  out.all.ls[[m]] <- cbind(lc.df, N_adult.df) %>%
    gather(year, N_adult, (ncol(lc.df)+1):ncol(.)) %>% 
    mutate(year=as.numeric(year),
           N_all=unlist(N_all.df),
           B=c(B.tot))
  out.df.ls[[m]] <- lc.df %>%
    mutate(N_adult.0=N_adult.tot[,1],
           N_all.0=N_all.tot[,1],
           B.0=B.tot[,1],
           N_adult.final=N_adult.tot[,dem_par$tmax+1],
           N_all.final=N_all.tot[,dem_par$tmax+1],
           B.final=B.tot[,dem_par$tmax+1])
}

# plots
library(viridis)
plot.dir <- "out/mgmt/"
gifs <- T
theme_set(theme_bw())
out.df <- dplyr::bind_rows(out.df.ls, .id="mgmt")
out.all <- dplyr::bind_rows(out.all.ls, .id="mgmt")
out.property <- out.all %>% filter(!is.na(Property)) %>%
  group_by(mgmt, year, Property) %>% 
  summarise(N_adult=sum(N_adult), N_all=sum(N_all), B=sum(B), nCell=n())
mgmt_comp.df <- out.df %>% 
  select(mgmt, lon, lat, id, id.in, in.UNH, Property, N_all.final) %>% 
  tidyr::spread(mgmt, N_all.final) %>%
  mutate(propDiff.none=(none-none)/none,
         propDiff.stated=(stated-none)/none,
         propDiff.reality=(reality-none)/none,
         propDiff.agg=(aggressive-none)/none) %>%
  select(-c(none, stated, reality, aggressive)) %>%
  gather(mgmt, prDiff, contains("propDiff"))
write.csv(out.df, paste0(plot.dir, res, "_out_df.csv"))
write.csv(out.all, paste0(plot.dir, res, "_out_all.csv"))
write.csv(out.property, paste0(plot.dir, res, "_out_property.csv"))
write.csv(mgmt_comp.df, paste0(plot.dir, res, "_mgmt_comp.csv"))
unh.bbox <- with(filter(lc.df, in.UNH), c(range(lon), range(lat)))

# abundance through time
ggplot(filter(out.all, in.UNH), aes(year, N_adult, group=id)) + 
  geom_line(alpha=0.5) + facet_grid(mgmt~Property)
ggplot(filter(out.all, in.UNH), aes(year, N_all, group=id)) + 
  geom_line(alpha=0.5) + facet_grid(mgmt~Property)
ggplot(filter(out.all, in.UNH), aes(year, B, group=id)) + 
  geom_line(alpha=0.5) + facet_grid(mgmt~Property)
ggplot(out.property, aes(year, N_adult/nCell, group=Property, colour=Property)) +
  geom_line() + facet_wrap(~mgmt) + 
  labs(y=paste("Mean buckthorn adult density per", res))
ggplot(out.property, aes(year, N_all/nCell, group=Property, colour=Property)) +
  geom_line() + facet_wrap(~mgmt) + 
  labs(y=paste("Mean buckthorn total density per", res))
ggplot(out.property, aes(year, B/nCell, group=Property, colour=Property)) +
  geom_line() + facet_wrap(~mgmt) + 
  labs(y=paste("Mean buckthorn seed bank density per", res))

# gifs
if(gifs) {
  library(gganimate)
  if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
  anim_save(paste0(plot.dir, "N_mgmt.gif"),
            animate(ggplot(filter(out.all, 
                                          lon>=unh.bbox[1] & lon<=unh.bbox[2] &
                                            lat>=unh.bbox[3] & lat<=unh.bbox[4]), 
                                   aes(lon, lat)) + 
                      geom_tile(aes(fill=log(N), colour=in.UNH), size=0.5) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(values=c(NA, "white")) +
                      transition_time(year) +  
                      facet_wrap(~mgmt) +
                      ggtitle("Adult log(abundance). Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
  anim_save(paste0(plot.dir, "B_mgmt.gif"),
            animate(ggplot(filter(out.all, 
                                  lon>=unh.bbox[1] & lon<=unh.bbox[2] &
                                    lat>=unh.bbox[3] & lat<=unh.bbox[4]), 
                           aes(lon, lat)) + 
                      geom_tile(aes(fill=log(B), colour=in.UNH), size=0.5) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(values=c(NA, "white")) +
                      transition_time(year) +  
                      facet_wrap(~mgmt) +
                      ggtitle("Seed bank log(abundance). Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
}



ggplot(mgmt_comp.df, aes(lon, lat, fill=prDiff)) + geom_tile() +
  scale_fill_gradient2(low="blue", high="red", limits=c(-1,1)) + facet_wrap(~mgmt)

ggplot(filter(out.df, 
              lon>=unh.bbox[1] & lon<=unh.bbox[2] &
                lat>=unh.bbox[3] & lat<=unh.bbox[4]),
       aes(lon, lat, fill=log(N_adult.final), colour=in.UNH)) +
  geom_tile() + scale_colour_manual(values=c(NA, "white")) +
  scale_fill_viridis(option="B") + facet_wrap(~mgmt)

ggplot(filter(mgmt_comp.df, 
              lon>=unh.bbox[1] & lon<=unh.bbox[2] &
                lat>=unh.bbox[3] & lat<=unh.bbox[4]), 
       aes(lon, lat, fill=prDiff, colour=in.UNH)) + 
  geom_tile() + scale_colour_manual(values=c(NA, "black")) +
  scale_fill_gradient2(low="blue", high="red", limits=c(-1,1)) + facet_wrap(~mgmt)

ggplot(filter(mgmt_comp.df, in.UNH)) + 
  geom_density(aes(x=prDiff, colour=mgmt)) +
  facet_wrap(~Property, scales="free_y")
ggplot(filter(mgmt_comp.df, in.UNH)) + 
  geom_boxplot(aes(x=mgmt, y=prDiff)) +
  facet_wrap(~Property)

