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
Packages <- c("here", "gbPopMod", "tidyverse", "magrittr", "viridis")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
set.seed(1)




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
#--- load libraries & landscape
gifs <- TRUE
res <- "9km2"  # 9km2 or 20ac
load(paste0("data/USDA_", res, ".rda"))  # loads dataframe as lc.df
parcel.df <- read.csv(paste0("data/USDA_", res, "_parcels.csv")) %>% 
  filter(inbd) %>% droplevels %>% mutate(id.pp=row_number())
ngrid <- nrow(lc.df)  # number of cells in bounding box
ncell <- sum(lc.df$inbd)  # number of inbound cells
npp <- nrow(parcel.df)  # number of pixel-parcels
LCs <- names(lc.df)[4:9]
pp.ls <- map(1:ncell, ~which(parcel.df$id.in==.))
coord.init <- c(739235.9, 4753487) # 1922: first record; herbarium_records.R
cell_side <- mean(diff(sort(unique(lc.df$lon))))
cell.init <- lc.df$id[which(abs(lc.df$lon-coord.init[1]) < cell_side/2 & 
                              abs(lc.df$lat-coord.init[2]) < cell_side/2)]
# Note that this works with three sets of indexes:
# - id: rectangular grid ID; identified for all cells
# - id.in: inbound ID; identified for inbound cells, NA if out of bounds
# - id.pp: pixel-parcel ID; identifier for each decision-making entity 
# This is to allow simpler spatial calculations (i.e., for SDD neighborhoods) 
# and for allocating buckthorn individuals and management effects appropriately



#--- set parameters
tmax <- 100
# global parameters: ?set_g_p
g.p <- set_g_p(tmax=tmax, n.cores=1,  
               K=c(5223180, 0, 770790, 770790, 770790, 770790),
               sdd.max=7, sdd.rate=1.4)



#--- FAKE ECONOMIC DECISION REGRESSIONS
get.actions2 <- function(N_i, K_i, PP_R, PP_T, TRT, N_EFF, T_EFF, 
                         W, CHEM, B, COST){ 
  
  WTP <- B[1] * W * PP_R + 
    B[2] * ((1-N_i*N_EFF/K_i)/(1-N_i/K_i)-1) * PP_R + 
    B[3] * CHEM * PP_R + 
    B[4] * N_i * T_EFF * PP_T
  
  return(TRT[which.max(WTP - COST)])
}

trt.i <- list(trt=c("N", "M", "C", "B"), # treatment name
              N_eff=c(0, 0.2, 0.3, 0.9), # effectiveness = mortality rate
              biomass_eff=c(0, 0.1, 0.125, 0.2), # effect on biomass
              cost=c(0, 200, 100, 300), # cost
              chem=c(0, 0, 1, 1), # chemical used?
              wild_eff=c(0, 1, 1, 1)) # effect on wildlife
WTP_b <- c(10, 40, -20, 60) # willingness-to-pay slopes
actions <- array(NA, dim=c(npp, tmax+1)) # store actions in each time step




#--- initialize SDD neighborhoods & buckthorn populations (just load eventually)
sdd <- sdd_set_probs(ncell, lc.df, g.p, verbose=T)
B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
N <- array(0, dim=c(npp, tmax+1, 6, max(g.p$m))) # [px-parcel, year, LC, age]
N.0.px <- pop_init(ngrid, g.p, lc.df, p.0=cell.init)[lc.df$inbd,,]
N[,1,,] <- round(N.0.px[parcel.df$id.in,,]*parcel.df$Grid_Proportion)
for(l in 1:6) { if(g.p$m[l] < 7) { N[,,l,g.p$m[l]:(max(g.p$m)-1)] <- NA } }
  
parcel.df$K_pp <- (as.matrix(parcel.df[,LCs]) %*% g.p$K) %>%
  multiply_by(parcel.df$Grid_Proportion) %>% c %>% round




system.time({  # just to get a sense
pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
for(k in 1:g.p$tmax) {
  ##---
  ## 2. Management plans are decided
  ##---
  # Economic decision model: which cells perform which management actions? 
  N.k <- rowSums(N[,k,,7])  # buckthorn adult abundance
  presence_pp <- which(N.k>0)  # id.pp with adult buckthorn
  action_pp <- sapply(presence_pp,
                      function(x) get.actions2(N_i=N.k[x],
                                               K_i=parcel.df$K_pp[x],
                                               PP_R=parcel.df$Wilderness_PP[x],
                                               PP_T=parcel.df$Timber_PP[x],
                                               TRT=trt.i$trt,
                                               N_EFF=trt.i$N_eff,
                                               T_EFF=trt.i$biomass_eff,
                                               W=trt.i$wild_eff,
                                               CHEM=trt.i$chem,
                                               B=WTP_b,
                                               COST=trt.i$cost))
  actions[presence_pp,k] <- action_pp # store pixel-parcel decisions
  mech_chem_id.pp <- presence_pp[action_pp != "N"] # id.pp that are treating
  
  # Set control parameters, identify treated pixel-parcels
  control_par <- set_control_p(null_ctrl=F,
                               man.i=mech_chem_id.pp, # mech/chem parcels
                               man.trt=setNames(trt.i$N_eff, trt.i$trt)) 
  # Specify which parcels use which treatments
  if(length(mech_chem_id.pp)>0) {
    mech_chem.i <- data.frame(id=mech_chem_id.pp,
                              Trt=action_pp[action_pp != "N"]) 
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
})




# plots
plot.dir <- "out/econ/"
gifs <- T
theme_set(theme_bw())
if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
N.tot <- t(sapply(1:ncell, function(x) apply(N[pp.ls[[x]],,,max(g.p$m)], 2, sum)))
N.df <- setNames(as.data.frame(N.tot), 1:ncol(N.tot))
N.df$id.in <- 1:nrow(N.df)
out.all <- left_join(lc.df, N.df, by="id.in") %>%
  gather(year, N, (ncol(lc.df)+1):ncol(.)) %>% 
  mutate(year=as.numeric(year),
         B=c(B))
actions.ls <- setNames(vector("list", 4), trt.i$trt)
for(i in 1:length(trt.i$trt)) {
  actions.ls[[i]] <- map_dfr(setNames(1:ncell, 1:ncell), 
                          ~colSums(actions[pp.ls[[.]],]==trt.i$trt[i], na.rm=T)/length(pp.ls[[.]])) %>%
    mutate(year=1:(g.p$tmax+1),
           action=trt.i$trt[i]) %>%
    gather(id.in, prop, 1:ncell) %>%
    mutate(id.in=as.numeric(id.in))
}
out.actions <- full_join(filter(lc.df, inbd), 
                         do.call("rbind", actions.ls), by="id.in")
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
                      geom_tile(aes(fill=log(N))) +
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
  anim_save(paste0(plot.dir, "Parcel_decisions.gif"),
            animate(ggplot(out.actions, aes(lon, lat)) + 
                      geom_tile(aes(fill=prop)) +
                      scale_fill_viridis(name="Proportion of\npixel-parcels", 
                                         option="B") + 
                      transition_time(year) +
                      facet_wrap(~action) +
                      ggtitle("Management decisions. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=1100, height=1000, units="px", fps=10))
}


