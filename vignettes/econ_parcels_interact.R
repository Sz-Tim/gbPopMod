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
  filter(inbd) %>% droplevels %>%
  mutate(id.pp=row_number())
ngrid <- nrow(lc.df)  # number of cells in bounding box
ncell <- sum(lc.df$inbd)  # number of inbound cells
npp <- nrow(parcel.df)  # number of pixel-parcels
LCs <- names(lc.df)[4:9]
pp.ls <- map(1:ncell, ~which(parcel.df$id.in==.))
# Note that this works with three sets of indexes:
# - id: rectangular grid ID; identified for all cells
# - id.in: inbound ID; identified for inbound cells, NA if out of bounds
# - id.pp: pixel-parcel ID; identifier for each decision-making entity 
# This is to allow simpler spatial calculations (i.e., for SDD neighborhoods) 
# and for allocating buckthorn individuals and management effects appropriately



#--- set parameters
tmax <- 150
# global parameters: ?set_g_p
g.p <- set_g_p(tmax=tmax, n.cores=1,  
               K=c(5223180, 0, 770790, 770790, 770790, 770790),
               sdd.max=7, sdd.rate=1.4)



#--- FAKE ECONOMIC DECISION REGRESSIONS
# pr(treat) = antilogit(b1 + b2*log(N_buckthorn_adults))
mech_chem_b <- c(-2, 0.5)
grd_cover_b <- c(-4, 0.6)
N_seq <- exp(seq(log(1), log(max(g.p$K)), length.out=100))
plot(log(N_seq), antilogit(cbind(1,log(N_seq)) %*% mech_chem_b), 
     ylim=c(0,1), type="l", ylab="Pr(manual treatment)", xlab="log(buckthorn)")
plot(log(N_seq), antilogit(cbind(1,log(N_seq)) %*% grd_cover_b), 
     ylim=c(0,1), type="l", ylab="Pr(ground cover)", xlab="log(buckthorn)")



#--- load or initialize SDD neighborhoods & buckthorn populations
# if(!file.exists("data/inits_sdd.rds")) {
  sdd <- sdd_set_probs(ncell, lc.df, g.p, verbose=T)
  # saveRDS(sdd, "data/inits_sdd.rds")
# } else {
#   sdd <- readRDS("data/inits_sdd.rds")
# }
# if(!file.exists("data/inits_B.rds")) {
  B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
  # saveRDS(B, "data/inits_B.rds")
# } else {
#   B <- readRDS("data/inits_B.rds")
# }
# if(!file.exists("data/inits_N.rds")) {
  N <- array(0, dim=c(npp, tmax+1, 6, max(g.p$m))) # [px-parcel, year, LC, age]
  N.0.px <- pop_init(ncell, g.p, lc.df, p.0=2000) # [cell, LC, age]
  N[,1,,] <- N.0.px[parcel.df$id.in,,]*parcel.df$Grid_Proportion
  for(l in 1:6) { if(g.p$m[l] < 7) { N[,,l,g.p$m[l]:(max(g.p$m)-1)] <- NA } }
  # saveRDS(N, "data/inits_N.rds")
# } else {
#   N <- readRDS("data/inits_N.rds")
# }
parcel.df$K_pp <- (as.matrix(parcel.df[,LCs]) %*% g.p$K) %>%
  multiply_by(parcel.df$Grid_Proportion) %>% round




system.time({
pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
for(k in 1:g.p$tmax) {
  ##---
  ## 2. Management plans are decided
  ##---
  # economic decision model
  #   which cells perform which management actions? Here, management is 
  #   probabilistic based on adult buckthorn abundance and the made up slopes
  #   shown above. The two treatment types are generated independently, though
  #   I believe in your model they'll be considered jointly.
  N.mx <- cbind(1, log(rowSums(N[,k,,7]) + 0.0001))
  mech_chem_id.pp <- which(rbinom(npp, 1, antilogit(N.mx %*% mech_chem_b))==1)
  grd_cover_id.pp <- which(rbinom(npp, 1, antilogit(N.mx %*% grd_cover_b))==1)
  
  # Set control parameters, identify treated pixel-parcels
  #   Here, there is one manual treatment (mechanical + chemical) and one
  #   ground treatment (cover crop = turf). Any additional treatments listed
  #   used in mech_chem.i or grd_cover.i just need to be established first via
  #   set_control_p(man.trt=c(MC=0.8, ...), grd.trt=c(Cov=0.005, ...))
  control_par <- set_control_p(null_ctrl=F,
                               man.i=mech_chem_id.pp, # mech/chem parcels
                               man.trt=c(MC=0.9, M=0.3), # mech/chem effectiveness
                               grd.i=grd_cover_id.pp, # ground cover parcels
                               grd.trt=c(Cov=0.005)) # ground cover effectiveness
  # Specify which parcels use which treatments
  if(length(mech_chem_id.pp)>0) {
    mech_chem.i <- data.frame(id=mech_chem_id.pp,
                              Trt=sample(c("MC", "M"), length(mech_chem_id.pp),
                                         replace=T)) # random mech vs mech+chem
  } else {mech_chem.i <- NULL}
  if(length(grd_cover_id.pp)>0) {
    grd_cover.i <- data.frame(id=grd_cover_id.pp,
                              Trt="Cov") # sets Trt="Cov" for all rows
  } else {grd_cover.i <- NULL}
  
  
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
                          sdd, control_par, grd_cover.i, mech_chem.i)
  N[,k+1,,] <- out$N
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
N.tot <- apply(N[,,,max(g.p$m)], 1:2, sum) # sum adults across LCs
N.tot <- t(sapply(1:ncell, function(x) apply(N.tot[pp.ls[[x]],], 2, sum)))
N.df <- setNames(as.data.frame(N.tot), 1:ncol(N.tot))
N.df$id.in <- 1:nrow(N.df)
out.all <- left_join(lc.df, N.df, by="id.in") %>%
  gather(year, N, (ncol(lc.df)+1):ncol(.)) %>% 
  mutate(year=as.numeric(year),
         B=c(B))
out.df <- lc.df %>%
  mutate(N.0=out.all$N[out.all$year==1],
         B.0=B[,1],
         N.final=out.all$N[out.all$year==g.p$tmax+1],
         B.final=B[,g.p$tmax+1])

# final maps
final.p <- ggplot(out.df, aes(lon, lat))
final.p + geom_tile(aes(fill=log(N.final))) + scale_fill_viridis(option="B")
final.p + geom_tile(aes(fill=log(B.final))) + scale_fill_viridis(option="B") 
final.p + geom_tile(aes(fill=N.final > 0))

# abundance through time
ggplot(out.all, aes(year, N, group=id)) + geom_line(alpha=0.5)

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
                    width=800, height=600, units="px", fps=5))
  anim_save(paste0(plot.dir, "Seed_abundance.gif"),
            animate(ggplot(out.all, aes(lon, lat)) + 
                      geom_tile(aes(fill=log(B))) +
                      scale_fill_viridis(option="B") + 
                      transition_time(year) +  
                      ggtitle("Seed abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px", fps=5))
}


