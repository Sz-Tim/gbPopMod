# This script runs simulations on the cluster
# The authorization key for the private repo is stored in the
# working directory in git_gb_token.R

# load most recent version of package on GitHub
# source("hpc/git_gb_token.R")
# devtools::install_github("Sz-Tim/gbPopMod", auth_token=git_gb_token, 
#                          dependencies=TRUE)

# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "stringr", "here", "doSNOW",
              "fastmatch", "scales", "gganimate")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

# set landscape file
lc_f <- "data/9km_car.csv"

# set parameters
g.p <- set_g_p(tmax=50, lc.r=100, lc.c=100, n.cores=4, 
               m=c(3,3,7,7,7,7), sdd.max=5, sdd.rate=1, N.p.t0=40)
par.ls <- set_sensitivity_pars(names(g.p)[10:25])
nSamp <- 300

# load landscape
lc.df <- read_csv(lc_f) %>% 
  filter(y >= (max(.$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.in=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)
N.init <- pop_init(ngrid, g.p, lc.df)

# run sensitivity analysis
out <- global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, lc.df, 
                          sdd.pr, N.init, control.p=NULL, verbose=T)

nMetric <- 8
nPar <- ncol(out$results)-nMetric
brt.sum <- vector("list", nMetric)
for(i in 1:nMetric) {
  metric <- names(out$results)[nPar+i]
  emulate_sensitivity(out$results, par.ls, g.p$n.cores, resp=metric)
  brt.sum[[i]] <- emulation_summary(metric)
}
ri.df <- map_dfr(brt.sum, ~.$ri.df)
cvDev.df <- map_dfr(brt.sum, ~.$cvDev.df)
betaDiv.df <- map_dfr(brt.sum, ~.$betaDiv.df)


# store output
write_csv(out$results, "out/sensitivity_results.csv")
write_csv(ri.df, "out/BRT_RI.csv")
write_csv(cvDev.df, "out/BRT_cvDev.csv")
write_csv(betaDiv.df, "out/BRT_betaDiv.csv")
#saveRDS(out, "out/sensitivity_out.rds")
