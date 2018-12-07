# This script runs simulations for the global sensitivity analysis. It is
# designed to run many simulations in parallel. The number of parallel
# simulations is controlled by assigning set_g_p(n.cores=) in line 23.

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=T,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")

## The GSA varies the following parameters simultaneously:
##  p.f: pr(flower)
##  mu: mean(fruits | flower)
##  gamma: mean(seeds/fruit)
##  m: maturation age
##  p.c: pr(fruit consumed by bird)
##  sdd.rate: 1/mean(short distance dispersal distance in cells)
##  sdd.max: maximum short distance dispersal distance
##  bird.hab: bird relative habitat preferences
##  n.ldd: long distance dispersal events per year
##  s.c: survival rate for seeds consumed by birds
##  s.B: survival rate in seed bank
##  s.M: survival rate for juveniles
##  s.N: survival rate for adults
##  K: carrying capacity
##  g.B: pr(germination from seed bank)
##  p: pr(establish | germination)
##  N.0: initial abundance in cell with first historical record

########
## Setup
########
# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW","fastmatch")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

# set parameters
res <- c("20ac", "9km2")[1]
g.p <- set_g_p(tmax=50, lc.r=Inf, lc.c=Inf, N.p.t0=1, n.cores=4)
par.ls <- set_sensitivity_pars(names(g.p)[10:26][-15], "gb", res)
par.ls$N.0 <- list(param="N.0", type="int", LC=0, min=1, max=100)
g.p$N.0 <- 10
nSamp <- 400

# load landscape
load(paste0("data/USDA_", res, ".rda")) # loads landscape as lc.df
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
cell.init <- get_pt_id(lc.df, c(739235.9, 4753487)) # 1922: herbarium_records.R


########
## Run model
########
# run sensitivity analysis
out.dir <- paste0("out/", res, "/")
global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, lc.df, 
                   sdd=NULL, cell.init, control.p=NULL, verbose=T, 
                   sim.dir=paste0(out.dir, "sims/"))

out <- dir(out.dir, recursive=T, full.names=T) %>% map_dfr(read.csv)



########
## Emulate output
########
nMetric <- 6  # pOcc, pSB, pK, meanNg0, medNg0, sdNg0
nPar <- ncol(out)-nMetric
brt.sum <- vector("list", nMetric)
for(i in 1:nMetric) {
  metric <- names(out)[nPar+i]
  emulate_sensitivity(out, par.ls, g.p$n.cores, resp=metric,
                      brt.dir=paste0(out.dir, "brt/"))
  brt.sum[[i]] <- emulation_summary(metric, paste0(out.dir, "brt/"))
}



########
## Store emulation results
########
write_csv(map_dfr(brt.sum, ~.$ri.df), paste0(out.dir, "BRT_RI.csv"))
write_csv(map_dfr(brt.sum, ~.$cvDev.df), paste0(out.dir, "BRT_cvDev.csv"))
write_csv(map_dfr(brt.sum, ~.$betaDiv.df), paste0(out.dir, "BRT_betaDiv.csv"))



