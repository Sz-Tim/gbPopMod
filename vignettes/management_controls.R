# Example of management treatments
# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "stringr", "here", "doSNOW",
              "fastmatch", "scales", "gganimate")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())
data(lc.rct)

# set parameters
n.sim <- 3
g.p <- set_g_p(tmax=150, lc.r=50, lc.c=50, n.cores=1)
control.p <- set_control_p(null_ctrl=FALSE, 
                           man.i=1300:1800,  # cells with manual controls
                           nTrt.man=NA,  # for random cell assignment
                           man.trt=c(M=0.05, C=0.3, MC=0.95),
                           grd.i=1300:1500, # cells with ground cover controls
                           nTrt.grd=NA,  # for random cell assignment
                           grd.trt=c(Lit=0.005, Cov=0.01, Com=0.00001),
                           chg.i=NULL  # cells with timber harvest
                           )

# land cover
lc.df <- lc.rct %>% 
  filter(y >= (max(lc.rct$y) - g.p$lc.r) & x <= g.p$lc.c) %>%
  mutate(id=row_number(), 
         id.inbd=min_rank(na_if(inbd*id, 0)))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# short distance dispersal neighborhoods
sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)

# initialize populations
N.init <- pop_init(ngrid, g.p, lc.df)

# run simulation
if(g.p$n.cores > 1) {
  require(doParallel)
  p.c <- makeCluster(g.p$n.cores)
  registerDoParallel(p.c)
  out.p <- foreach(i=1:n.sim) %dopar% {
    gbPopMod::run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, 
                      N.init, control.p, verbose=F)
  }
  stopCluster(p.c)
  out.ad <- map(out.p, ~.$N[,,max(g.p$age.f)]) %>% unlist %>% 
    array(., dim=c(ngrid, g.p$tmax+1, n.sim))
  out.sb <- map(out.p, ~.$N.sb) %>% unlist %>% 
    array(., dim=c(ngrid, g.p$tmax+1, n.sim))
} else {
  out.ad <- array(dim=c(ngrid, g.p$tmax+1, n.sim))
  out.sb <- array(dim=c(ngrid, g.p$tmax+1, n.sim))
  for(s in 1:n.sim) {
    out <- run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, 
                   control.p, verbose=T)
    out.ad[,,s] <- out$N[,,max(g.p$age.f)]
    out.sb[,,s] <- out$N.sb
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
p.wd <- paste0("out/", ncell, "_t", g.p$tmax, "_", g.p$edges, "/")
if(!dir.exists(here::here("out"))) dir.create(here::here("out"))
if(!dir.exists(here::here(p.wd))) dir.create(here::here(p.wd))

# save plots
save_pars(p.wd, g.p, control.p)
make_plots_final_t(p.wd, g.p, filter(N.out, year==max(N.out$year)))
make_plots_gifs(p.wd, g.p, N.out)
make_plots_lc(p.wd, lc.df)



