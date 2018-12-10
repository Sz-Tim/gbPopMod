# This script sets up the landscape and initial populations. The first record in
# the study region is in 1922, and so one iteration is run for 96 years
# (1922-2018) with the single initial population. The average buckthorn
# abundance in each cell in 2018 is then used as the initial state for
# simulations beginning in the present day.

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=F,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")

########
## Setup
########
library(gbPopMod); library(doSNOW); library(foreach)

# set parameters
n_sim <- 100
res <- c("20ac", "9km2")[1]
dem_par <- set_g_p(tmax=96, n.cores=10)
if(res == "9km2") {
  dem_par$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  dem_par$sdd.max <- 7
  dem_par$sdd.rate <- 1.1
  dem_par$n.ldd <- 3
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


# store SDD neighborhoods
if(!dir.exists("data/inits/temp")) dir.create("data/inits/temp", recursive=T)
saveRDS(sdd, paste0("data/inits/sdd_", res, ".rds"))


N_0 <- pop_init(ngrid, dem_par, lc.df, p.0=cell.init, N.0=10)
N.tmax <- array(0, dim=c(ngrid, 6, max(dem_par$m), n_sim))
B.tmax <- matrix(0, nrow=ngrid, ncol=n_sim)

p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
sim.out <- foreach(s=1:n_sim,
                   .packages=c("gbPopMod", "stringr")) %dopar% {
  N.k <- N_0
  B.k <- rep(0, ngrid)
  for(k in 1:(dem_par$tmax-1)) {
    out <- iterate_pop(ngrid, ncell, N.k, B.k, dem_par, lc.df, sdd)
    N.k <- out$N
    B.k <- out$B
  }
  s.f <- str_pad(s, 4, "left", "0")
  saveRDS(N.k, paste0("data/inits/temp/", res, "_N_", s.f, ".rds"))
  saveRDS(B.k, paste0("data/inits/temp/", res, "_B_", s.f, ".rds"))
  return(s)
}
stopCluster(p.c)



# store abundances
N.tmax <- dir("data/inits/temp", paste0(res, "_N_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
B.tmax <- dir("data/inits/temp", paste0(res, "_B_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
saveRDS(N.tmax, paste0("data/inits/N_2018_", res, ".rds"))
saveRDS(B.tmax, paste0("data/inits/B_2018_", res, ".rds"))


# calculate extra SDD info
p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
sdd.ji.rows <- foreach(x=1:ncell) %dopar% { which(sdd$sp.df$j.idin==x) }
stopCluster(p.c)
sdd.ji <- parallel::mclapply(sdd.ji.rows, function(x) sdd$sp.df$i.idin[x]) 
p.ji <- parallel::mclapply(sdd.ji.rows, function(x) sdd$sp.df$pr[x]) 
saveRDS(sdd.ji, paste0("data/inits/sdd_ji_", res, ".rds"))
saveRDS(p.ji, paste0("data/inits/p_ji_", res, ".rds"))















