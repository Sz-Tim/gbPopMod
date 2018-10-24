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
tmax <- 50

#--- load landscape, initial populations
Packages <- c("here", "gbPopMod", "tidyverse", "magrittr", "viridis")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
load(paste0("data/USDA_", res, ".rda"))  # loads dataframe as lc.df
lc.df$Rgrw <- 0
lc.df <- lc.df[,c(1:4,15,5:14)]
parcel.df <- read.csv(paste0("data/USDA_", res, "_parcels.csv")) %>% 
  filter(inbd) %>% droplevels %>% mutate(id.pp=row_number(), Rgrw=0)
parcel.df <- parcel.df[c(1:11,23,12:22)]
N_init <- readRDS(paste0("data/inits/N_2018_", res, ".rds"))
N_0 <- array(0, dim=dim(N_init)+c(0,1,0))
N_0[,c(1,3:7),] <- N_init
rm(N_init)
B_0 <- readRDS(paste0("data/inits/B_2018_", res, ".rds"))
sdd <- readRDS(paste0("data/inits/sdd_", res, ".rds"))
sdd.ji <- readRDS(paste0("data/inits/sdd_ji_", res, ".rds"))
ngrid <- nrow(lc.df)  # number of cells in bounding box
ncell <- sum(lc.df$inbd)  # number of inbound cells
npp <- nrow(parcel.df)  # number of pixel-parcels
LCs <- names(lc.df)[4:10]
pp.ls <- map(1:ncell, ~which(parcel.df$id.in==.))



#--- set parameters
# global demographic parameters: ?set_g_p
g.p <- set_g_p(tmax=tmax, n.cores=1)  
if(res=="9km2") {
  g.p$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  g.p$sdd.max <- 7
  g.p$sdd.rate <- 1.4
}
# adjust parameters for Regrow category: same parameters as Open
# Categories: 1=Opn, 2=Rgrw, 3=Oth, 4=Dec, 5=Evg, 6=WP, 7=Mxd
g.p$n.lc <- 7
lc.par <- which(map_lgl(g.p, ~length(.)>1))
g.p[lc.par] <- map(g.p[lc.par], ~.[c(1,1:6)])
for(l in 1:g.p$n.lc) {
  parcel.df[[paste0("K_", LCs[l])]] <- ceiling(parcel.df[[LCs[l]]] * g.p$K[l] *
                                                 parcel.df$Grid_Proportion)
}
parcel.df$K_pp <- rowSums(select(parcel.df, contains("K_")))
parcel.df_orig <- parcel.df  # store original proportions if forest regrowth 



#--- FAKE ECONOMIC & TIMBER HARVEST REGRESSIONS
get.actions2 <- function(N_i, K_i, PP_R, PP_T, TRT, N_EFF, T_EFF, 
                         W, CHEM, B, COST) { 
  WTP <- B[1] * W * PP_R + 
    B[2] * ((1-N_i*N_EFF/K_i)/(1-N_i/K_i)-1) * PP_R + 
    B[3] * CHEM * PP_R + 
    B[4] * N_i * T_EFF * PP_T
  return(TRT[which.max(WTP - COST)])
}
time_to_forest <- function(T.0, beta, K, N.t0, N.t1) {
  (T.0 - 1) * (1 + beta * N.t1/K)/(1 + beta * N.t0/K)
}



#--- FAKE FOREST HARVEST DECISION PARAMETERS
nharvest <- 1000  # number of pixel-parcels that harvest in a given year
beta <- 1
L.0 <- 10
forest.pp <- which(names(parcel.df) %in% c("Dec", "Evg", "WP", "Mxd"))
cut.df <- data.frame(id.pp=numeric(0), T_regrow=numeric(0))



#--- FAKE MANAGEMENT DECISION PARAMETERS
trt.i <- list(trt=c("N", "M", "C", "B"), # treatment name
              N_eff=c(0, 0.3, 0.5, 0.9), # effectiveness = mortality rate
              biomass_eff=c(0, 0.1, 0.15, 0.3), # effect on biomass
              cost=c(0, 100, 200, 500), # cost
              chem=c(1, 1, -1, -1), # chemical used? (-1 = YES, 1 = No)
              wild_eff=c(0, 1, 1, 1)) # effect on wildlife
WTP_b <- c(10, 40, -20, 60) # willingness-to-pay slopes
Mgmt1_LC <- LCs[c(1,4:7)]  # Opn, Dec, Evg, WP, Mxd
Mgmt2_LC <- LCs[2]  # Rgrw
act1.mx <- act2.mx <- array("N", dim=c(npp, tmax)) # store actions in each year



#--- initialize SDD neighborhoods & buckthorn populations (just load eventually)
B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
B[,1] <- B_0
N <- array(0, dim=c(npp, tmax+1, 7, max(g.p$m))) # [px-parcel, year, LC, age]
N.0.px <- N_0[lc.df$inbd,,]
N[,1,,] <- round(N.0.px[parcel.df$id.in,,]*parcel.df$Grid_Proportion)
N[,1,,7] <- pmin(N[,1,,7], as.matrix(parcel.df$K_pp * parcel.df[,LCs]))
  




########################################
## SIMULATION LOOP
####################

pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
for(k in 1:g.p$tmax) {
  
  ##---
  ## 2. Change land cover: Timber harvest & regrowth on managed properties
  ##---
  # Forest harvest model: which pixel-parcels harvest forest?
  ## NOTE: You will need to specify how much of each forest type is cleared
  ##    For regrowth, you'll need to calculate & track time-until-forest for 
  ##    each harvested pixel-parcel. 
  
  # harvest: id.pp, proportion cleared for each forest type (assuming 100% here)
  # NOTE: K and N are both just for within landcover category 'Rgrw': N[,,2,]
  if(k>1 && nrow(cut.df)>0) {
    cut.df$T_regrow <- time_to_forest(T.0=cut.df$T_regrow,
                                      beta=beta,
                                      K=parcel.df$K_Rgrw[cut.df$id.pp],
                                      N.t0=N[cut.df$id.pp,k-1,2,7],
                                      N.t1=N[cut.df$id.pp,k,2,7])
    cut_pool <- (1:npp)[-cut.df$id.pp]
  } else {
    cut_pool <- 1:npp
  }
  cut_k <- data.frame(id.pp=sample(cut_pool, nharvest),
                      Dec=rep(1, nharvest),
                      Evg=rep(1, nharvest),
                      WP=rep(1, nharvest),
                      Mxd=rep(1, nharvest))
  # adjust N: forest becomes Rgrw
  N_cut <- N[cut_k$id.pp,k,4:7,]
  for(i in 1:7) { # iterate across ages; multiply N by % cut
    N_cut[,,i] <- N_cut[,,i] * as.matrix(cut_k[,-1])
  }
  N[cut_k$id.pp,k,2,] <- apply(N_cut, c(1,3), sum)
  N[cut_k$id.pp,k,4:7,] <- N[cut_k$id.pp,k,4:7,] - N_cut
  # regrowth
  if(any(cut.df$T_regrow <= 0)) {
    regrow_id.pp <- cut.df$id.pp[cut.df$T_regrow <= 0]
    cut.df <- filter(cut.df, T_regrow > 0)
    # regrow to original proportions
    parcel.df[regrow_id.pp, LCs] <- parcel.df_orig[regrow_id.pp, LCs]
    regrow_id <- sort(unique(parcel.df$id[regrow_id.pp]))
    # adjust N: Rgrw becomes forest
    N_regrow <- N[regrow_id.pp,k,2,]
    for(i in 4:7) {
      N[regrow_id.pp,k,i,] <- N_regrow * parcel.df[regrow_id.pp, LCs[i]]
    }
    N[regrow_id.pp,k,2,] <- 0
  } else {
    regrow_id.pp <- regrow_id <- NULL
  }
  altered_id.pp <- unique(c(cut_k$id.pp, regrow_id.pp))
  # convert cut forest percentages to Rgrw
  parcel.df[cut_k$id.pp, "Rgrw"] <- rowSums(parcel.df[cut_k$id.pp, forest.pp] * 
                                              cut_k[,-1])
  # reduce forest cover by specified proportions
  parcel.df[cut_k$id.pp, forest.pp] <- parcel.df[cut_k$id.pp, forest.pp] *
    (1 - cut_k[,-1])
  # udpate lc.df
  cut_id <- sort(unique(parcel.df$id[cut_k$id.pp]))
  altered_id <- sort(unique(c(regrow_id, cut_id)))
  altered_LC <- parcel.df %>% filter(id %in% altered_id) %>%
    group_by(id) %>%
    summarise_at(LCs, sum)
  altered_LC[,-1] <- altered_LC[,-1]/rowSums(altered_LC[,-1])
  lc.df[altered_id, LCs] <- altered_LC[,-1]
  # udpate SDD neighborhoods
  sdd.alter <- tibble(id.in=unique(unlist(sdd.ji[lc.df$id.in[altered_id]])),
                      id=lc.df$id[match(id.in, lc.df$id.in)])
  sdd_new <- sdd_update_probs(lc.df, g.p, sdd.alter, sdd$i, lc.col=4:10)
  sdd$i[,,1,sdd.alter$id.in] <- sdd_new$i
  sdd$sp[sdd.alter$id.in] <- sdd_new$sp
  # update carrying capacities
  for(l in 1:g.p$n.lc) {
    parcel.df[altered_id.pp, paste0("K_", LCs[l])] <- ceiling(
      parcel.df[altered_id.pp,LCs[l]] * g.p$K[l] * 
        parcel.df$Grid_Proportion[altered_id.pp])
  }
  parcel.df$K_pp[altered_id.pp] <- rowSums(
    parcel.df[altered_id.pp, paste0("K_", LCs)])
  # append current year harvests to cut.df
  cut.df <- bind_rows(cut.df,
                      data.frame(id.pp=cut_k$id.pp,
                                 T_regrow=(1+beta*N[cut_k$id.pp,k,2,7]/
                                             parcel.df$K_Rgrw[cut_k$id.pp])*L.0))
  
  
  ##---
  ## 3. Management plans are decided
  ##---
  # Economic decision model: which pixel-parcels perform which mgmt actions? 
  # N[id.pp (1:npp), year (1:tmax), LC (1:7), age (1:7)]
  N.k.1 <- rowSums(N[,k,which(LCs %in% Mgmt1_LC),7])
  N.k.2 <- N[,k,which(LCs %in% Mgmt2_LC),7]
  presence_pp.1 <- which(N.k.1>0)  # id.pp with adult buckthorn
  presence_pp.2 <- which(N.k.2>0)  # id.pp with adult buckthorn
  action_pp.1 <- sapply(presence_pp.1,
                        function(x) get.actions2(N_i=N.k.1[x],
                                                 K_i=sum(parcel.df[x,paste0("K_", Mgmt1_LC)]),
                                                 PP_R=parcel.df$Wilderness_PP[x],
                                                 PP_T=parcel.df$Timber_PP[x],
                                                 TRT=trt.i$trt,
                                                 N_EFF=trt.i$N_eff,
                                                 T_EFF=0,
                                                 W=trt.i$wild_eff,
                                                 CHEM=trt.i$chem,
                                                 B=WTP_b,
                                                 COST=trt.i$cost))
  action_pp.2 <- sapply(presence_pp.2,
                        function(x) get.actions2(N_i=N.k.2[x],
                                                 K_i=sum(parcel.df[x,paste0("K_", Mgmt2_LC)]),
                                                 PP_R=parcel.df$Wilderness_PP[x],
                                                 PP_T=parcel.df$Timber_PP[x],
                                                 TRT=trt.i$trt,
                                                 N_EFF=trt.i$N_eff,
                                                 T_EFF=trt.i$biomass_eff,
                                                 W=trt.i$wild_eff,
                                                 CHEM=trt.i$chem,
                                                 B=WTP_b,
                                                 COST=trt.i$cost))
  act1.mx[presence_pp.1,k] <- action_pp.1 # store pixel-parcel decisions
  act2.mx[presence_pp.2,k] <- action_pp.2 
  mech_chem_id.pp <- unique(c(presence_pp.1[action_pp.1 != "N"], 
                              presence_pp.2[action_pp.2 != "N"]))

  # Set control parameters, identify treated pixel-parcels
  control_par <- set_control_p(null_ctrl=F,
                               man.trt=setNames(trt.i$N_eff, trt.i$trt))
  # Specify which parcels use which treatments
  if(length(mech_chem_id.pp)>0) {
    # If treatment varies by land cover type, this needs to be a dataframe with
    # nLC + 1 columns (id, Trt.Opn, Trt.Oth, etc) and a row per pixel-parcel. 
    mech_chem.i <- data.frame(id=mech_chem_id.pp,
                              Trt.Opn=act1.mx[mech_chem_id.pp,k],
                              Trt.Rgrw=act2.mx[mech_chem_id.pp,k],
                              Trt.Oth="N",
                              Trt.Dec=act1.mx[mech_chem_id.pp,k],
                              Trt.Evg=act1.mx[mech_chem_id.pp,k],
                              Trt.WP=act1.mx[mech_chem_id.pp,k],
                              Trt.Mxd=act1.mx[mech_chem_id.pp,k])
  } else {
    mech_chem.i <- NULL
  }
  
  
  ##---
  ## 4. MANAGEMENT IS IMPLEMENTED, BUCKTHORN GROWS & SPREADS
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
  anim_save(paste0(plot.dir, "Parcel_decisions_1.gif"),
            animate(ggplot(out.actions, aes(lon, lat)) + 
                      geom_tile(aes(fill=prop.1)) +
                      scale_fill_viridis(name="Proportion of\npixel-parcels", 
                                         option="B") + 
                      transition_time(year) +
                      facet_wrap(~action) +
                      ggtitle(paste("Management decisions for",
                                    Mgmt1_LC, "Year {frame_time}")),
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
                                    Mgmt2_LC, "Year {frame_time}")),
                    nframes=n_distinct(out.all$year), 
                    width=1100, height=1000, units="px", fps=10))
}


