# Example of management treatments
# load libraries
Packages <- c("here", "doSNOW", "fastmatch", "scales", "gganimate",
              "gbPopMod", "tidyverse", "magrittr")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

theme_set(theme_bw())

# set parameters
n.sim <- 4
g.p <- set_g_p(tmax=50, n.cores=4, sdd.max=5, sdd.rate=1.3, N.p.t0=1)
control.p <- set_control_p(null_ctrl=FALSE, 
                           t.trt=20,
                           man.i=1300:1800,  # cells with manual controls
                           pTrt.man=NA,  # for random cell assignment
                           man.trt=c(M=0.05, C=0.3, MC=0.95),
                           grd.i=1300:1500, # cells with ground cover controls
                           pTrt.grd=NA,  # for random cell assignment
                           grd.trt=c(Lit=0.005, Cov=0.01, Com=0.00001),
                           lc.chg=FALSE,
                           pChg=.1,
                           chg.i=NULL  # cells with timber harvest
                           )

# land cover
lc.df <- read_csv("data/USDA_9km2.csv") # USDA_9km2.csv or USDA_20ac.csv
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)
id.i <- lc.df %>% select(id, id.in)

# short distance dispersal neighborhoods
sdd.pr <- sdd_set_probs(ncell, lc.df, g.p, verbose=T)

# initialize populations
N.init <- pop_init(ngrid, g.p, lc.df)


# how is this related to run_sim? This simulation does not have management controls
# out.lam <- run_sim_lambda(ngrid, ncell, g.p, lambda_values, 
#                          lc.df, sdd.pr$i, N.init, TRUE)
# TODO: is this needed?
# out <- lc.df %>% mutate(lam=c(as.matrix(lc.df[,4:9]) %*% lambda_values),
#                        N=out.lam$N[,g.p$tmax+1])

# run simulation
if(g.p$n.cores > 1) {
  require(doSNOW)
  p.c <- makeCluster(g.p$n.cores)
  registerDoSNOW(p.c)
  out.p <- foreach(i=1:n.sim, .packages="gbPopMod") %dopar% {
    run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, control.p, verbose=F)
  }
  stopCluster(p.c)
  out.ad <- map(out.p, ~.$N[,,max(g.p$m)]) %>% unlist %>% 
    array(., dim=c(ngrid, g.p$tmax+1, n.sim))
  out.sb <- map(out.p, ~.$B) %>% unlist %>% 
    array(., dim=c(ngrid, g.p$tmax+1, n.sim))
} else {
  out.ad <- out.sb <- array(dim=c(ngrid, g.p$tmax+1, n.sim))
  for(s in 1:n.sim) {
    out <- run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, control.p)
    out.ad[,,s] <- out$N[,,max(g.p$m)]
    out.sb[,,s] <- out$B
    rm(out)
    cat("Finished simulation", s, "of", n.sim, "\n")
  }
}

# plots
ad.mn <- apply(out.ad, 1:2, mean)
N.out <- cbind(lc.df, ad.mn) %>% as.tibble %>%
  gather(year, ad_Ab, (1:ncol(ad.mn)) + ncol(lc.df)) %>%
  mutate(ad_sd=c(apply(out.ad, 1:2, sd)),
         ad_pP=c(apply(out.ad > 0, 1:2, mean)),
         sb_Ab=c(apply(log(out.sb+1), 1:2, mean)),
         sb_sd=c(apply(log(out.sb+1), 1:2, sd)),
         sb_pP=c(apply(out.sb > 0, 1:2, mean)),
         ad_L5=c(ad.mn > 0))
N.out$ad_L5 <- N.out$ad_L5 + c(ad.mn > 5)
N.out$year <- str_pad(N.out$year, 3, "left", "0")

# prepare plot paths & filenames
p.wd <- paste0("out/", ncell, "_t", g.p$tmax, "/")
if(!dir.exists(here::here(p.wd))) dir.create(here::here(p.wd), recursive=T)

# save plots
save_pars(p.wd, g.p, control.p)
make_plots_final_t(p.wd, g.p, filter(N.out, year==max(N.out$year)))
make_plots_gifs(p.wd, g.p, N.out)
make_plots_lc(p.wd, lc.df)



