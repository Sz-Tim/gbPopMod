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
tmp.dir <- "data/inits/temp/"
if(!dir.exists(tmp.dir)) dir.create(tmp.dir, recursive=T)
saveRDS(sdd, paste0("data/inits/sdd_", res, ".rds"))


N_0 <- pop_init(ngrid, dem_par, lc.df, p.0=cell.init, N.0=10)

p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
sim.out <- foreach(s=1:n_sim,
                   .packages=c("gbPopMod", "stringr")) %dopar% {
  out <- run_sim(ngrid, ncell, dem_par, lc.df, sdd, N_0, NULL, F, dem_par$tmax)
  s.f <- str_pad(s, 4, "left", "0")
  saveRDS(out$N[,1,], paste0(tmp.dir, res, "_N_", s.f, ".rds"))
  saveRDS(out$B[,1], paste0(tmp.dir, res, "_B_", s.f, ".rds"))
  saveRDS(out$nSd[,1], paste0(tmp.dir, res, "_nSd_", s.f, ".rds"))
  saveRDS(out$nSdStay[,1], paste0(tmp.dir, res, "_nSdStay_", s.f, ".rds"))
  saveRDS(out$D[,1], paste0(tmp.dir, res, "_D_", s.f, ".rds"))
  return(s)
}
stopCluster(p.c)



# store summarized output
N.tmax <- dir(tmp.dir, paste0(res, "_N_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
B.tmax <- dir(tmp.dir, paste0(res, "_B_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
nSd.tmax <- dir(tmp.dir, paste0(res, "_nSd_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
nSdStay.tmax <- dir(tmp.dir, paste0(res, "_nSdStay_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
D.tmax <- dir(tmp.dir, paste0(res, "_D_"), full.names=T) %>%
  map(readRDS) %>%
  Reduce(`+`, .)/n_sim
saveRDS(N.tmax, paste0("data/inits/N_2018_", res, ".rds"))
saveRDS(B.tmax, paste0("data/inits/B_2018_", res, ".rds"))
saveRDS(nSd.tmax, paste0("data/inits/nSd_2018_", res, ".rds"))
saveRDS(nSdStay.tmax, paste0("data/inits/nSdStay_2018_", res, ".rds"))
saveRDS(D.tmax, paste0("data/inits/D_2018_", res, ".rds"))


# calculate extra SDD info
p.c <- makeCluster(dem_par$n.cores); registerDoSNOW(p.c)
sdd.ji.rows <- foreach(x=1:ncell) %dopar% { which(sdd$sp.df$j.idin==x) }
stopCluster(p.c)
sdd.ji <- parallel::mclapply(sdd.ji.rows, function(x) sdd$sp.df$i.idin[x]) 
p.ji <- parallel::mclapply(sdd.ji.rows, function(x) sdd$sp.df$pr[x]) 
saveRDS(sdd.ji, paste0("data/inits/sdd_ji_", res, ".rds"))
saveRDS(p.ji, paste0("data/inits/p_ji_", res, ".rds"))















