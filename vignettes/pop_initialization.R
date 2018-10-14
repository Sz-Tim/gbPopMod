# This script sets up the landscape and initial populations. The first record in
# the study region is in 1922, and so one iteration is run for 96 years
# (1922-2018) with the single initial population. The average buckthorn
# abundance in each cell in 2018 is then used as the initial state for
# simulations beginning in the present day.

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=T,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")

########
## Setup
########
library(gbPopMod); library(doSNOW); library(foreach)

# set parameters
n_sim <- 3
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
# sdd <- list(sp=sdd_set_probs(ncell, lc.df, dem_par)$sp)
sdd <- sdd_set_probs(ncell, lc.df, dem_par)
p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
sdd.ji.rows <- foreach(x=1:ncell) %dopar% { which(sdd$sp.df$j.idin==x) }
stopCluster(p.c)
sdd.ji <- lapply(sdd.ji.rows, function(x) sdd$sp.df$i.idin[x]) 
p.ji <- lapply(sdd.ji.rows, function(x) sdd$sp.df$pr[x]) 
N_0 <- pop_init(ngrid, dem_par, lc.df, p.0=cell.init, N.0=10)

N.tmax <- array(0, dim=c(ngrid, 6, max(dem_par$m), n_sim))
B.tmax <- matrix(0, nrow=ngrid, ncol=n_sim)

pb <- txtProgressBar(min=0, max=n_sim*dem_par$tmax, width=80, style=3)
for(s in 1:n_sim) {
  N.k <- N_0
  B.k <- rep(0, ngrid)
  for(k in 1:(dem_par$tmax-1)) {
    out <- iterate_pop(ngrid, ncell, N.k, B.k, dem_par, lc.df, sdd)
    N.k <- out$N
    B.k <- out$B
    setTxtProgressBar(pb, (s-1)*dem_par$tmax + k)
  }
  N.tmax[,,,s] <- N.k
  B.tmax[,s] <- B.k
}
close(pb)


# store initial abundances & SDD neighborhoods
if(!dir.exists("data/inits/")) dir.create("data/inits/", recursive=T)
saveRDS(apply(N.tmax, 1:3, mean), paste0("data/inits/N_2018_", res, ".rds"))
saveRDS(rowMeans(B.tmax), paste0("data/inits/B_2018_", res, ".rds"))
saveRDS(sdd, paste0("data/inits/sdd_", res, ".rds"))
saveRDS(sdd.ji, paste0("data/inits/sdd_ji_", res, ".rds"))
saveRDS(p.ji, paste0("data/inits/p_ji_", res, ".rds"))














