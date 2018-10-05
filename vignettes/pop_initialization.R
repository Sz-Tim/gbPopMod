# This script sets up the landscape and initial populations. The first record in
# the study region is in 1922, and so one iteration is run for 96 years
# (1922-2018) with the single initial population. The buckthorn distribution in
# 2018 is then used as the initial state for simulations beginning in the
# present day.

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
res <- c("20ac", "9km2")[2]
dem_par <- set_g_p(tmax=96)
if(res == "9km2") {
  dem_par$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  dem_par$sdd.max <- 7
  dem_par$sdd.rate <- 1.4
}

# load landscape
load(paste0("data/USDA_", res, ".rda")) # loads landscape as lc.df
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# 1922: first record; herbarium_records.R
coord.init <- c(739235.9, 4753487)
cell.init <- get_pt_id(lc.df, coord.init)




########
## Simulate buckthorn
########
sdd <- list(sp=sdd_set_probs(ncell, lc.df, dem_par)$sp)
N_0 <- pop_init(ngrid, dem_par, lc.df, p.0=cell.init)

# run to year 96
N <- array(0, dim=c(ngrid, dem_par$tmax+1, 6, max(dem_par$m)))
N[,1,,] <- N_0

# B = seed bank; dim=[cell, year]
B <- matrix(0, nrow=ngrid, ncol=dem_par$tmax+1)

for(k in 1:dem_par$tmax) {
  out <- iterate_pop(ngrid, ncell, N[,k,,], B[,k], dem_par, lc.df, sdd, 
                     NULL, NULL, NULL)
  
  # update abundances
  N[,k+1,,] <- out$N
  B[,k+1] <- out$B
}

# store initial abundances & SDD neighborhoods
saveRDS(N[,96,,], "data/N_2018.rds")
saveRDS(B[,96], "data/B_2018.rds")
saveRDS(sdd, "data/sdd.rds")














