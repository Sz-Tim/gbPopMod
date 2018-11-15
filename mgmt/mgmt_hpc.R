# This script runs simulations of different management plans for the University
# of New Hampshire woodlands, which include 22 properties totaling nearly 3,500
# acres. The code is written to run iterations in parallel on a high performance
# cluster, and the number of cores to use can be specified with
# set_g_p(n.cores=...). We explore four different broad plans:
# 1) No management
# 2) Current management plan (property-specific)
# 3) Manual removal for populations  > 500 individuals/20-acres
# 4) Manual removal + cover crop for populations > 500 individuals/20-acres

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
dem_par <- set_g_p(tmax=20)
if(res == "9km2") {
  dem_par$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  dem_par$sdd.max <- 7
  dem_par$sdd.rate <- 1.4
}

# load landscape
# load(paste0("data/USDA_", res, ".rda")) # loads landscape as lc.df
lc.df <- read.csv("data/USDA_UNH_mgmt.csv")
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initial populations
N_0 <- readRDS(paste0("data/inits/N_2018_", res, ".rds"))[lc.df$id.full,,]
B_0 <- readRDS(paste0("data/inits/B_2018_", res, ".rds"))[lc.df$id.full]
sdd <- sdd_set_probs(ncell, lc.df, dem_par, verbose=T)

# management plans
mgmt_plans <- dir("mgmt", "plan", full.names=T) %>% map(read.csv) %>%
  setNames(str_remove(str_remove(dir("mgmt", "plan"), "plan_"), ".csv"))
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
    mgmt[[m]] <- list(ctrl_par=set_control_p(null_ctrl=FALSE,
                                             man.trt=c("M"=0.4, "C"=0.9, 
                                                       "MC"=0.96, "N"=0)),
                      manual.i=mgmt_m[,grep("id|Trt", names(mgmt_m))],
                      thresh=mgmt_m$thresh,
                      trt.int=mgmt_m$trt.cycle,
                      trt.followup=mgmt_m$trt.followup)
  }
} 
out.df.ls <- out.all.ls <- setNames(vector("list", 4), names(mgmt))

for(m in seq_along(mgmt)) {
  # storage objects
  # N = abundances; dim=[cell, year, LC, age]
  N <- array(0, dim=c(ngrid, dem_par$tmax+1, 6, max(dem_par$m)))
  N[,1,,] <- N_0
  # dem_par$m = age at maturity; varies by LC, so some are NA
  # N[,,, dem_par$m] stores the number of buckthorn adults
  # N[,,, 1:(dem_par$m-1)] stores the number of buckthorn juveniles 
  
  # B = seed bank; dim=[cell, year]
  B <- matrix(0, nrow=ngrid, ncol=dem_par$tmax+1)
  B[,1] <- B_0
  
  # ids of all managed cells
  managed <- unique(c(mgmt[[m]]$manual.i$id, mgmt[[m]]$cover.i$id))
  
  # iterate through years
  for(k in 1:dem_par$tmax) {
    
    # decide which cells are managed this year
    if(names(mgmt)[m] != "none" && 
       any((k %% mgmt[[m]]$trt.int) <= mgmt[[m]]$trt.followup)) {  
      managed_k <- managed[apply(N[managed,k,,], 1, sum) > mgmt[[m]]$thresh &
                             (k %% mgmt[[m]]$trt.int) <= mgmt[[m]]$trt.followup]
      manual_k <- mgmt[[m]]$manual.i[mgmt[[m]]$manual.i$id %in% managed_k,]
      cover_k <- mgmt[[m]]$cover.i[mgmt[[m]]$cover.i$id %in% managed_k,]
    } else {
      manual_k <- NULL
      cover_k <- NULL
    }
    
    # buckthorn management, population growth, dispersal
    out <- iterate_pop(ngrid, ncell, N[,k,,], B[,k], dem_par, lc.df[,-(15:18)], 
                       sdd, mgmt[[m]]$ctrl_par, cover_k, manual_k)
    
    # update abundances
    N[,k+1,,] <- out$N
    B[,k+1] <- out$B
    k <- k+1
  }
  
  # summarise
  N.tot <- apply(N[,,,max(dem_par$m)], 1:2, sum) # sum adults across LCs
  N.df <- as.data.frame(N.tot); names(N.df) <- 1:ncol(N.df)
  out.all.ls[[m]] <- cbind(lc.df, N.df) %>%
    gather(year, N, (ncol(lc.df)+1):ncol(.)) %>% 
    mutate(year=as.numeric(year),
           B=c(B))
  out.df.ls[[m]] <- lc.df %>%
    mutate(N.0=N.tot[,1],
           B.0=B[,1],
           N.final=N.tot[,dem_par$tmax+1],
           B.final=B[,dem_par$tmax+1])
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
  summarise(N=sum(N), B=sum(B), nCell=n())
# abundance through time
ggplot(filter(out.all, in.UNH), aes(year, log(N), group=id, colour=Property)) + 
  geom_line() + facet_wrap(~mgmt)
ggplot(filter(out.all, in.UNH), aes(year, log(B), group=id, colour=Property)) + 
  geom_line() + facet_wrap(~mgmt)
ggplot(out.property, aes(year, N/nCell, group=Property, colour=Property)) +
  geom_line() + facet_wrap(~mgmt) + 
  labs(y=paste("Buckthorn adult density per", res))
ggplot(out.property, aes(year, B/nCell, group=Property, colour=Property)) +
  geom_line() + facet_wrap(~mgmt) + 
  labs(y=paste("Buckthorn seed bank density per", res))

# gifs
if(gifs) {
  library(gganimate)
  if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
  anim_save(paste0(plot.dir, "N_mgmt.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=log(N), colour=in.UNH), size=0.5) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(values=c(NA, "white")) +
                      transition_time(year) +  
                      facet_wrap(~mgmt) +
                      ggtitle("Adult log(abundance). Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
  anim_save(paste0(plot.dir, "B_mgmt.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=log(B), colour=in.UNH), size=0.5) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(values=c(NA, "white")) +
                      transition_time(year) +  
                      facet_wrap(~mgmt) +
                      ggtitle("Seed bank log(abundance). Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
}



