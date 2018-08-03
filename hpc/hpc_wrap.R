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

# set resolution & parameter span
res <- c("20ac", "9km2")[2]
par_span <- c("total", "gb")[2]

# set parameters
g.p <- set_g_p(tmax=30, lc.r=Inf, lc.c=Inf, n.cores=4, N.p.t0=4)
par.ls <- set_sensitivity_pars(names(g.p)[10:26][-(6:7)], par_span, res)
nSamp <- 400

# load landscape
lc.df <- read_csv(paste0("data/USDA_", res, ".csv")) %>% 
  filter(y >= (max(.$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.in=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
N.init <- pop_init(ngrid, g.p, lc.df)

# run sensitivity analysis
out <- global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, lc.df, 
                          sdd=NULL, N.init, control.p=NULL, verbose=T, 
                          sim.dir=paste0("out/", res, "/", par_span, "/sims/"))
write_csv(out, paste0("out/", res, "/", par_span, "/gsa_results.csv"))

nMetric <- 8
nPar <- ncol(out)-nMetric
brt.sum <- vector("list", nMetric)
for(i in 1:nMetric) {
  metric <- names(out)[nPar+i]
  emulate_sensitivity(out, par.ls, g.p$n.cores, resp=metric, 
                      brt.dir=paste0("out/", res, "/", par_span, "/brt/"))
  brt.sum[[i]] <- emulation_summary(metric, 
                                    paste0("out/", res, "/", par_span, "/brt/"))
}

write_csv(map_dfr(brt.sum, ~.$ri.df), 
          paste0("out/", res, "/", par_span, "/BRT_RI.csv"))
write_csv(map_dfr(brt.sum, ~.$cvDev.df), 
          paste0("out/", res, "/", par_span, "/BRT_cvDev.csv"))
write_csv(map_dfr(brt.sum, ~.$betaDiv.df), 
          paste0("out/", res, "/", par_span, "/BRT_betaDiv.csv"))



