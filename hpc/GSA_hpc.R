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
g.p <- set_g_p(tmax=50, lc.r=Inf, lc.c=Inf, n.cores=8, N.p.t0=1)
par.ls <- set_sensitivity_pars(names(g.p)[10:26][-9], par_span, res)
par.ls$N.0 <- list(param="N.0", type="int", LC=1, min=rep(1,6), max=par.ls$K$max)
g.p$N.0 <- round(g.p$K/2)
nSamp <- 100

# load landscape
load(paste0("data/USDA_", res, ".rda"))
lc.df <- lc.df %>% 
  filter(y >= (max(.$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.in=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
coord.init <- c(739235.9, 4753487) # 1922: first record; herbarium_records.R
cell_side <- mean(diff(sort(unique(lc.df$lon))))
cell.init <- lc.df$id[which(abs(lc.df$lon-coord.init[1]) < cell_side/2 & 
                              abs(lc.df$lat-coord.init[2]) < cell_side/2)]

# run sensitivity analysis
out <- global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, lc.df, 
                          sdd=NULL, cell.init, control.p=NULL, verbose=T, 
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



