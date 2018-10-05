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
Packages <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW","fastmatch")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

# set parameters
dem_par <- set_g_p(tmax=50, lc.r=Inf, lc.c=Inf, n.cores=4, 
                   sdd.max=7, sdd.rate=1.4)

# load landscape
load(paste0("data/USDA_9km2.rda")) # loads landscape as lc.df
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# load UNH woodland properties
library(measurements); library(sf)
unh.mgmt <- readxl::read_xlsx("~/Desktop/UNH_woodlands_summary.xlsx", 1) %>%
  mutate(Lat=gsub("°", " ", str_remove(Lat, "'")),
         Lon=gsub("°", " ", str_remove(Lon, "'"))) %>%
  mutate(Lat=as.numeric(conv_unit(Lat, from='deg_dec_min', to='dec_deg')),
         Lon=as.numeric(conv_unit(Lon, from='deg_dec_min', to='dec_deg'))) %>%
  st_as_sf(coords=c("Lon", "Lat"), crs=4326) %>%
  st_transform(crs=32618)
unh.coord <- st_coordinates(unh.mgmt)
unh.id <- sapply(1:nrow(unh.coord), function(x) get_pt_id(lc.df, unh.coord[x,]))

# initial populations
N_0 <- readRDS("data/mgmt_N_2018.rds")
B_0 <- readRDS("data/mgmt_B_2018.rds")
sdd <- readRDS("data/mgmt_sdd.rds")

# management plans
mgmt <- list(none=list(ctrl_par=set_control_p(null_ctrl=TRUE),
                       manual.i=NULL,
                       cover.i=NULL),
             current=list(ctrl_par=set_control_p(null_ctrl=FALSE,
                                                 man.trt=c("MC"=0.9)),
                          manual.i=data.frame(id=unh.id,
                                              Trt="MC"),
                          cover.i=NULL,
                          thresh=1000,
                          trt.int=3),
             manual=list(ctrl_par=set_control_p(null_ctrl=FALSE,
                                                man.trt=c("C"=0.6)),
                         manual.i=data.frame(id=unh.id,
                                             Trt="C"),
                         cover.i=NULL,
                         thresh=500,
                         trt.int=2),
             manual_cover=list(ctrl_par=set_control_p(null_ctrl=FALSE,
                                                      man.trt=c("MC"=0.9),
                                                      grd.trt=c("Cov"=0.001)),
                               manual.i=data.frame(id=unh.id,
                                                   Trt="MC"),
                               cover.i=data.frame(id=unh.id,
                                                  Trt="Cov"),
                               thresh=500,
                               trt.int=5))
out.df.ls <- out.all.ls <- setNames(vector("list", 4), names(mgmt))

for(m in seq_along(mgmt)) {
  # storage objects
  # N = abundances; dim=[cell, year, LC, age]
  N <- array(0, dim=c(ngrid, dem_par$tmax+1, 6, max(dem_par$m)))
  N[,1,,] <- N_0
  # dem_par$m = age at maturity; varies by LC, so some are NA
  # N[,,, dem_par$m] stores the number of buckthorn adults
  # N[,,, 1:(dem_par$m-1)] stores the number of buckthorn juveniles
  for(l in 1:6) { 
    if(dem_par$m[l] < 7) N[,,l,dem_par$m[l]:(max(dem_par$m)-1)] <- NA
  }
  
  # B = seed bank; dim=[cell, year]
  B <- matrix(0, nrow=ngrid, ncol=dem_par$tmax+1)
  B[,1] <- B_0
  
  # ids of all managed cells
  managed <- unique(c(mgmt[[m]]$manual.i$id, mgmt[[m]]$cover.i$id))
  
  # iterate through years
  for(k in 1:dem_par$tmax) {
    
    # decide which cells are managed this year
    if(m != 1 && (k %% mgmt[[m]]$trt.int) == 0) {  
      managed_k <- managed[apply(N[managed,k,,], 1, sum) > mgmt[[m]]$thresh]
      manual_k <- mgmt[[m]]$manual.i[mgmt[[m]]$manual.i$id %in% managed_k,]
      cover_k <- mgmt[[m]]$cover.i[mgmt[[m]]$cover.i$id %in% managed_k,]
    } else {
      manual_k <- NULL
      cover_k <- NULL
    }
    
    # buckthorn management, population growth, dispersal
    out <- iterate_pop(ngrid, ncell, N[,k,,], B[,k], dem_par, lc.df, sdd, 
                       mgmt[[m]]$ctrl_par, cover_k, manual_k)
    
    # update abundances
    N[,k+1,,] <- out$N
    B[,k+1] <- out$B
  }
  # plots
  plot.dir <- "out/mgmt/"
  gifs <- T
  theme_set(theme_bw())
  if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
  N.tot <- apply(N[,,,max(dem_par$m)], 1:2, sum) # sum adults across LCs
  N.df <- as.data.frame(N.tot); names(N.df) <- 1:ncol(N.df)
  out.all.ls[[m]] <- cbind(lc.df, N.df) %>%
    gather(year, N, (ncol(lc.df)+1):ncol(.)) %>% 
    mutate(year=as.numeric(year)) %>%
    mutate(managed=id %in% unh.id)
  out.df.ls[[m]] <- lc.df %>%
    mutate(N.0=N.tot[,1],
           B.0=B[,1],
           N.final=N.tot[,dem_par$tmax+1],
           B.final=B[,dem_par$tmax+1],
           managed=id %in% unh.id)
}

out.df <- dplyr::bind_rows(out.df.ls, .id="mgmt")
out.all <- dplyr::bind_rows(out.all.ls, .id="mgmt")
# final maps
library(viridis)
final.p <- ggplot(out.df, aes(lon, lat)) + facet_wrap(~mgmt)
# final.p + geom_tile(aes(fill=log(N.final), colour=managed), size=0.5) + 
#   scale_fill_viridis(option="B") + scale_colour_manual(values=c(NA, "blue"))
# final.p + geom_tile(aes(fill=log(B.final), colour=managed), size=0.5) + 
#   scale_fill_viridis(option="B") + scale_colour_manual(values=c(NA, "blue"))

# abundance through time
ggplot(filter(out.all, managed), aes(year, N, group=id)) + 
  geom_line(alpha=0.5) + facet_wrap(~mgmt)

# gifs
if(gifs) {
  library(gganimate)
  anim_save(paste0(plot.dir, "N_mgmt.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=N, colour=managed), size=0.5) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(values=c(NA, "white")) +
                      transition_time(year) +  
                      facet_wrap(~mgmt) +
                      ggtitle("Adult abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
}



