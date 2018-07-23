# This script runs simulations on the cluster
# The authorization key for the private repo is stored in the
# working directory in git_gb_token.R

# load most recent version of package on GitHub
# source("hpc/git_gb_token.R")
# devtools::install_github("Sz-Tim/gbPopMod", auth_token=git_gb_token)

# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "stringr", "here", "doSNOW",
              "fastmatch", "scales")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

# set landscape file
lc_f <- "data/9km_car.csv"

# set parameters
g.p <- set_g_p(tmax=50, lc.r=100, lc.c=100, n.cores=16, 
               m=c(3,3,7,7,7,7), sdd.max=5, sdd.rate=1, N.p.t0=40)
par.ls <- set_sensitivity_pars(names(g.p)[10:25])
nSamp <- 10000

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
out <- global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, lc.df, sdd.pr, 
                          N.init, control.p=NULL, verbose=T, 
                          sim.dir=paste0("out/", par_span, "/sims/"))
write_csv(out, paste0("out/", par_span, "/sensitivity_results.csv"))

nMetric <- 8
nPar <- ncol(out)-nMetric
brt.sum <- vector("list", nMetric)
for(i in 1:nMetric) {
  metric <- names(out)[nPar+i]
  emulate_sensitivity(out, par.ls, g.p$n.cores, resp=metric, 
                      brt.dir=paste0("out/", par_span, "/brt/"))
  brt.sum[[i]] <- emulation_summary(metric, paste0("out/", par_span, "/brt/"))
}

write_csv(map_dfr(brt.sum, ~.$ri.df), 
          paste0("out/", par_span, "/BRT_RI.csv"))
write_csv(map_dfr(brt.sum, ~.$cvDev.df), 
          paste0("out/", par_span, "/BRT_cvDev.csv"))
write_csv(map_dfr(brt.sum, ~.$betaDiv.df), 
          paste0("out/", par_span, "/BRT_betaDiv.csv"))



