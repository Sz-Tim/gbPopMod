# This script runs simulations on the cluster
# The authorization key for the private repo is stored in the
# working directory in git_gb_token.R

# load most recent version of package on GitHub
# source("hpc/git_gb_token.R")
# devtools::install_github("Sz-Tim/gbPopMod", auth_token=git_gb_token)

# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "stringr", "here", "doSNOW",
              "fastmatch", "scales", "gganimate")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

# load files
par_i <- read.csv("data/param_ranges.csv", stringsAsFactors=F)
data(lc.rct)

# set parameters
g.p <- set_g_p(tmax=30, lc.r=200, lc.c=200, n.cores=3, sdd.max=10, N.p.t0=200)
pars <- par_i$p[c(7,11,12,16:18)]
nSamp <- 12

# munge
pars.rng <- filter(par_i, p %in% pars)
pars.rng <- pars.rng[match(pars, pars.rng$p), ]
lc.df <- lc.rct %>% 
  filter(y >= (max(lc.rct$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.in=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)
N.init <- pop_init(ngrid, g.p, lc.df)

# run sensitivity analysis
out <- global_sensitivity(pars, pars.rng, nSamp, ngrid, ncell, g.p, lc.df, 
                          sdd.pr, N.init, control.p=NULL, verbose=T)

# store output
saveRDS(out, "out/sensitivity_out.rds")
