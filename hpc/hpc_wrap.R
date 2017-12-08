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
n.sim <- 6
g.p <- set_g_p(tmax=10, lc.r=30, lc.c=30, n.cores=3)
control.p <- set_control_p()
p <- readRDS("hpc/p.rds")
p.seq <- readRDS("hpc/p_seq.rds")

# land cover
lc.df <- lc.rct %>% 
  filter(y >= (max(lc.rct$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.inbd=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# sensitivity loop
map2(p, p.seq, 
     ~run_sensitivity(.x, .y, n.sim, ngrid, ncell, g.p, control.p, lc.df))

