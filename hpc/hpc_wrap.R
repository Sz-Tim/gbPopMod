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
data(lc.rct)


# set parameters
n.sim <- 3
g.p <- set_g_p(tmax=5, lc.r=30, lc.c=30, n.ldd=1, n.cores=3)
control.p <- set_control_p()
p <- c("pr.sb", "pr.s.bird", "pr.s", "pr.f", "fec")
p.seq <- list(seq(0.1, 0.8, length.out=3),
              seq(0.3, 0.8, length.out=3),
              expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              expand_cnpy(Op=c(0.2, 0.9), Cl=c(0.2, 0.9), length_out=2),
              expand_cnpy(Op=c(100, 200), Cl=c(10, 40), length_out=2))

# land cover
lc.df <- lc.rct %>% 
  filter(y >= (max(lc.rct$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.inbd=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# population
set.seed(225)
N.init <- pop_init(ngrid, g.p, lc.df)

##---
## run sensitivity loop
##---

map2(p, p.seq, ~run_sensitivity(.x, .y, n.sim, ngrid, ncell, g.p, control.p,
                                lc.df, N.init))





