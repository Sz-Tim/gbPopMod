# USDA-NIFA Buckthorn model
# Tim Szewczyk

# This is an example of how the buckthorn model can interact with the economic
# decision model where decisions are made (and management actions taken) by each
# unique pixel-parcel. Here, 'pixel' refers to a grid cell on the land scape. A
# pixel may be owned by multiple parcels. 

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=T,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")
# The following packages are called by various gbPopMod functions:
# - here: easier file directory navigation
# - doSNOW: running in parallel
# - foreach: running in parallel
# - fastmatch: faster SDD neighborhood calculation
# - tidyverse: tidyr, dplyr, ggplot, purrr, tibble, stringr, forcats, readr
# - magrittr: additional pipe (%>%) functions
# - *scales: generating color scales in some plotting functions (PLOTS ONLY)
# - *gganimate: generating gifs (PLOTS ONLY)
# - *viridis: for pretty colors (PLOTS ONLY)
# - *dismo: boosted regression trees for sensitivity analysis (SENSITIVITY ONLY)

# Here, the economic model is written in R for illustration and the buckthorn
# population is iterated each time step. The basic flow is:
# 1. Set up and initialization. Store initial buckthorn population, etc
# 2. Management plans are decided
# 3. Management is implemented, buckthorn grows and spreads
# 4. Repeat 2-3 for tmax years
#
# Buckthorn dynamics occur within each pixel of the landscape, but pixels are
# potentially subdivided into parcels. Each pixel-parcel decides whether and 
# how to manage buckthorn on their property. These pixel-parcels are nonspatial,
# but control a specified proportion of the pixel, and we assume that the land
# cover proportions are constant within each pixel such that the pixel-parcels
# within a given pixel contain the same land cover composition.
#
# NOTE: Pixel = grid cell = cell
#
# This simulation works with three sets of indexes:
# - id: rectangular grid ID; identified for all cells
# - id.in: inbound ID; identified for inbound cells, NA if out of bounds
# - id.pp: pixel-parcel ID; identifier for each decision-making entity 
#
# We accommodate this difference in spatial resolution as follows:
# - K[id.pp]: Buckthorn carrying capacity per pixel-parcel = K[id] * % of cell
# - N[id.pp]: Track buckthorn abundance per pixel-parcel, constrained by K
# - Fruits[id]: Aggregate total fruit production per pixel
# - Seeds[id]: Disperse total seeds per pixel, uniform placement in target cell
# - Seed bank[id]: Total seed bank per pixel
# - Germinants[id]: Total germinating seeds per pixel
# - Seedlings[id.pp]: pixel-parcel = Germinants[id] * p_estab[id.pp] * % of cell
# - Management[id.pp]: based on N[id.pp], affect N[id.pp]



##---
## 1. SET UP AND INITIALIZATION
##---
#--- set general parameters
gifs <- TRUE
res <- "9km2"  # 9km2 or 20ac
tmax <- 10

#--- load landscape, initial populations
Packages <- c("here", "gbPopMod", "tidyverse", "magrittr", "viridis")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
load(paste0("data/USDA_", res, ".rda"))  # loads dataframe as lc.df
parcel.df <- read.csv(paste0("data/USDA_", res, "_parcels.csv")) %>% 
  filter(inbd) %>% droplevels %>% mutate(id.pp=row_number())
mgmt.lookup <- readxl::read_xlsx("econ/Full Combined Sheet with Decision Rules.xlsx", sheet=1) %>% 
  filter(!(LB==1 & UB==1)) %>% rename(mortality=`Removal Rate`)
open.lookup <- filter(mgmt.lookup, LC=="Open")
forest.lookup <- filter(mgmt.lookup, LC=="Forest")
N_init <- readRDS(paste0("data/inits/N_2018_", res, ".rds"))
N_0 <- array(0, dim=dim(N_init))
N_0 <- N_init
rm(N_init)
B_0 <- readRDS(paste0("data/inits/B_2018_", res, ".rds"))
sdd <- readRDS(paste0("data/inits/sdd_", res, ".rds"))
sdd.ji <- readRDS(paste0("data/inits/sdd_ji_", res, ".rds"))
ngrid <- nrow(lc.df)  # number of cells in bounding box
ncell <- sum(lc.df$inbd)  # number of inbound cells
npp <- nrow(parcel.df)  # number of pixel-parcels
LCs <- names(lc.df)[4:9]
pp.ls <- map(1:ncell, ~which(parcel.df$id.in==.))



#--- set parameters
# global demographic parameters: ?set_g_p
g.p <- set_g_p(tmax=tmax, n.cores=4)  
if(res=="9km2") {
  g.p$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  g.p$sdd.max <- 7
  g.p$sdd.rate <- 1.4
}
for(l in 1:g.p$n.lc) {
  parcel.df[[paste0("K_", LCs[l])]] <- ceiling(parcel.df[[LCs[l]]] * g.p$K[l] *
                                                 parcel.df$Grid_Proportion)
}
parcel.df$K_pp <- rowSums(select(parcel.df, contains("K_")))



#--- set management LCs
mgmt.ls <- vector("list", tmax)
Mgmt1_LC <- LCs[3:6]  # Dec, Evg, WP, Mxd
Mgmt2_LC <- LCs[1]  # Opn
lookup_mortality <- function(lookup.df, # lookup table for LC 
                             parcel.i, # single row of parcel.df
                             NK, # N/K for the parcel
                             Subsidy, # subsidy
                             k) {  # year
  # if it's a harvest year for this pixel
  if(parcel.i$MES==1 && k %% parcel.i$Harv_Freq) {
    mort_i <- filter(lookup.df, MES==parcel.i$MES & 
                       NMES==parcel.i$NMES &
                       S==Subsidy & 
                       NK > LB & 
                       NK <= UB)
    mort_harvest <- filter(mort_i, Harvest_Year==1)
    mort_nonharvest <- filter(mort_i, Harvest_Year==0)
    
    if(nrow(mort_harvest)==0) mort_harvest[1,] <- 0
    if(nrow(mort_nonharvest)==0) mort_nonharvest[1,] <- 0
    trt_i <- (1-parcel.i$Perc_Harv) * mort_nonharvest[1,8:11] +
      parcel.i$Perc_Harv * mort_harvest[1,8:11]
    
  } else {
    mort_i <- filter(lookup.df, MES==parcel.i$MES & 
                       NMES==parcel.i$NMES &
                       Harvest_Year==0 &
                       S==Subsidy & 
                       NK > LB & 
                       NK <= UB)
    if(nrow(mort_i)==0) mort_i[1,] <- 0
    trt_i <- mort_i[1,8:11]
  }
  
  return(unlist(trt_i))
}



#--- initialize SDD neighborhoods & buckthorn populations (just load eventually)
B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
B[,1] <- B_0
N <- array(0, dim=c(npp, tmax+1, g.p$n.lc, max(g.p$m))) # [px-parcel, year, LC, age]
N.0.px <- N_0[lc.df$inbd,,]
N[,1,,] <- round(N.0.px[parcel.df$id.in,,]*parcel.df$Grid_Proportion)
N[,1,,7] <- pmin(N[,1,,7], as.matrix(parcel.df$K_pp * parcel.df[,LCs]))
  




########################################
## SIMULATION LOOP
####################

pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
for(k in 1:g.p$tmax) {
  
  ##---
  ## 2. Management plans are decided
  ##---
  # Economic decision model: which pixel-parcels perform which mgmt actions? 
  # N[id.pp (1:npp), year (1:tmax), LC (1:7), age (1:7)]
  # calculate N/K for Forests (*.1) and for Open (*.2)
  NK.k.1 <- round(rowSums(N[,k,3:6,7])/ # Forests
                    rowSums(parcel.df[paste0("K_", LCs[3:6])]), 2)
  NK.k.2 <- round(N[,k,1,7]/parcel.df$K_Opn, 2) # Open
  presence_pp.1 <- which(NK.k.1>0)  # id.pp with adult buckthorn
  presence_pp.2 <- which(NK.k.2>0)  # id.pp with adult buckthorn
  
  # lookup treatment mortalities
  mort_pp.1 <- map(presence_pp.1, 
                   ~lookup_mortality(forest.lookup, parcel.df[.,], 
                                     NK.k.1[.], Subsidy=0.5, k))
  mort_pp.1 <- do.call('rbind', mort_pp.1)
  mort_pp.2 <- map(presence_pp.2, 
                   ~lookup_mortality(open.lookup, parcel.df[.,], 
                                     NK.k.2[.], Subsidy=0.5, k))
  mort_pp.2 <- do.call('rbind', mort_pp.2)
  
  # Set control parameters, identify treated pixel-parcels
  control_par <- set_control_p(null_ctrl=F)
  # Specify which parcels use which treatments
  mgmt.ls[[k]] <- data.frame(year=k,
                             id=1:npp,
                             NK.Forest=NK.k.1,
                             NK.Open=NK.k.2,
                             mort.Forest=0,
                             mort.Open=0,
                             use_mech.Forest=0,
                             use_chem.Forest=0,
                             use_both.Forest=0,
                             use_mech.Open=0,
                             use_chem.Open=0,
                             use_both.Open=0)
  # Store mortality rates and treatment methods
  if(any(c(mort_pp.1[,1], mort_pp.2[,1])>0)) {
    mort_non0.1 <- which(mort_pp.1[,1]>0)
    mort_non0.2 <- which(mort_pp.2[,1]>0)
    mech_chem.i <- list(id.trt.forest=data.frame(id=presence_pp.1[mort_non0.1],
                                                 mort=mort_pp.1[mort_non0.1, 1]),
                        id.trt.open=data.frame(id=presence_pp.2[mort_non0.2],
                                               mort=mort_pp.2[mort_non0.2, 1]))
    mgmt.ls[[k]]$mort.Forest[presence_pp.1[mort_non0.1]] <- mort_pp.1[mort_non0.1,1]
    mgmt.ls[[k]]$use_mech.Forest[presence_pp.1[mort_non0.1]] <- mort_pp.1[mort_non0.1,2]
    mgmt.ls[[k]]$use_chem.Forest[presence_pp.1[mort_non0.1]] <- mort_pp.1[mort_non0.1,3]
    mgmt.ls[[k]]$use_both.Forest[presence_pp.1[mort_non0.1]] <- mort_pp.1[mort_non0.1,4]
    mgmt.ls[[k]]$mort.Open[presence_pp.2[mort_non0.2]] <- mort_pp.2[mort_non0.2,1]
    mgmt.ls[[k]]$use_mech.Open[presence_pp.2[mort_non0.2]] <- mort_pp.2[mort_non0.2,2]
    mgmt.ls[[k]]$use_chem.Open[presence_pp.2[mort_non0.2]] <- mort_pp.2[mort_non0.2,3]
    mgmt.ls[[k]]$use_both.Open[presence_pp.2[mort_non0.2]] <- mort_pp.2[mort_non0.2,4]
  } else {
    mech_chem.i <- NULL
  }
  
  
  ##---
  ## 3. MANAGEMENT IS IMPLEMENTED, BUCKTHORN GROWS & SPREADS
  ##---
  # buckthorn population model
  # Parameters are updated for cells that performed management actions, then
  # buckthorn grows and spreads. This is implemented as:
  # 1. Update parameters based on management actions
  # 2. Adults produce fruit
  # 3. Fruit is consumed and dispersed, or dropped
  # 4. Seeds germinate and establish
  g.p$n.ldd <- ifelse(k %% 2 == 0, 1, 0)
  out <- iterate_pop_econ(parcel.df, pp.ls, N[,k,,], B[,k], g.p, lc.df, 
                          sdd, control_par, grd_cover.i=NULL, mech_chem.i)
  N.1g0 <- which(rowSums(out$N[,1,])>0)
  N[N.1g0,k+1,,] <- out$N[N.1g0,,]
  B[,k+1] <- out$B
  
  setTxtProgressBar(pb, k)
}
close(pb)




# plots
plot.dir <- "econ/"
gifs <- T
theme_set(theme_bw())
if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
N.tot <- t(sapply(1:ncell, function(x) apply(N[pp.ls[[x]],,,7], 2, sum)))
N.df <- setNames(as.data.frame(N.tot), 1:ncol(N.tot))
N.df$id.in <- 1:nrow(N.df)
out.all <- left_join(lc.df, N.df, by="id.in") %>%
  gather(year, N, (ncol(lc.df)+1):ncol(.)) %>% 
  mutate(year=as.numeric(year),
         B=c(B))
mgmt.all <- do.call('rbind', mgmt.ls)
act1.ls <- act2.ls <- setNames(vector("list", 4), trt.i$trt)
for(i in 1:length(trt.i$trt)) {
  act1.ls[[i]] <- map_dfr(setNames(1:ncell, 1:ncell), 
                          ~colSums(act1.mx[pp.ls[[.]],]==trt.i$trt[i])/
                           length(pp.ls[[.]])) %>%
    mutate(year=1:(g.p$tmax),
           action=trt.i$trt[i]) %>%
    gather(id.in, prop.1, 1:ncell) %>%
    mutate(id.in=as.numeric(id.in))
  act2.ls[[i]] <- map_dfr(setNames(1:ncell, 1:ncell), 
                          ~colSums(act2.mx[pp.ls[[.]],]==trt.i$trt[i])/
                            length(pp.ls[[.]])) %>%
    mutate(year=1:(g.p$tmax),
           action=trt.i$trt[i]) %>%
    gather(id.in, prop.2, 1:ncell) %>%
    mutate(id.in=as.numeric(id.in))
}
out.actions <- full_join(filter(lc.df, inbd), 
                         do.call("rbind", act1.ls), by="id.in") 
out.actions$prop.2 <- full_join(filter(lc.df, inbd), 
                                do.call("rbind", act2.ls), by="id.in")$prop.2
out.df <- lc.df %>%
  mutate(N.0=out.all$N[out.all$year==1],
         B.0=B[,1],
         N.final=out.all$N[out.all$year==g.p$tmax+1],
         B.final=B[,g.p$tmax+1])

# final abundance maps
final.p <- ggplot(out.df, aes(lon, lat))
final.p + geom_tile(aes(fill=log(N.final))) + scale_fill_viridis(option="B")
final.p + geom_tile(aes(fill=log(B.final))) + scale_fill_viridis(option="B") 

# abundance through time
ggplot(out.all, aes(year, log(N+1), group=id)) + geom_line(alpha=0.5)

# gifs
if(gifs) {
  library(gganimate)
  anim_save(paste0(plot.dir, "Adult_abundance.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=N)) +
                      scale_fill_viridis(option="B") + 
                      transition_time(year) +  
                      ggtitle("Adult abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px", fps=10))
  anim_save(paste0(plot.dir, "Seed_abundance.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=log(B))) +
                      scale_fill_viridis(option="B") + 
                      transition_time(year) +  
                      ggtitle("Seed abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px", fps=10))
  anim_save(paste0(plot.dir, "Parcel_decisions_1.gif"),
            animate(ggplot(out.actions, aes(lon, lat)) + 
                      geom_tile(aes(fill=prop.1)) +
                      scale_fill_viridis(name="Proportion of\npixel-parcels", 
                                         option="B") + 
                      transition_time(year) +
                      facet_wrap(~action) +
                      ggtitle(paste("Management decisions for",
                                    paste(Mgmt1_LC, collapse=" "), 
                                    "Year {frame_time}")),
                    nframes=n_distinct(out.all$year), 
                    width=1100, height=1000, units="px", fps=10))
  anim_save(paste0(plot.dir, "Parcel_decisions_2.gif"),
            animate(ggplot(out.actions, aes(lon, lat)) + 
                      geom_tile(aes(fill=prop.2)) +
                      scale_fill_viridis(name="Proportion of\npixel-parcels", 
                                         option="B") + 
                      transition_time(year) +
                      facet_wrap(~action) +
                      ggtitle(paste("Management decisions for",
                                    paste(Mgmt1_LC, collapse=" "),
                                    "Year {frame_time}")),
                    nframes=n_distinct(out.all$year), 
                    width=1100, height=1000, units="px", fps=10))
}


