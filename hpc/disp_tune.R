# This script runs simulations for the dispersal parameter tuning. It is
# designed to run many simulations in parallel. The number of parallel
# simulations is controlled by assigning set_g_p(n.cores=) in line 23.

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=F,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")

## The dispersal parameter tuning varies the following parameters simultaneously:
##  sdd.rate: 1/mean(short distance dispersal distance in cells)
##  sdd.max: maximum short distance dispersal distance
##  n.ldd: long distance dispersal events per year


########
## Setup
########
# load libraries
Packages <- c("gbPopMod", "tidyverse", "magrittr", "here", "foreach","doSNOW")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))

# set parameters
plots <- FALSE
res <- c("20ac", "9km2")[1]
g.p <- set_g_p(tmax=96, lc.r=Inf, lc.c=Inf, N.p.t0=1, n.cores=8)
par.ls <- set_sensitivity_pars(names(g.p)[c(15,16,18)], "gb", res)
g.p$N.0 <- 10
nSamp <- 2500
if(res=="9km2") {
  g.p$K <- c(3133908, 0, 462474, 462474, 462474, 462474)
  g.p$sdd.max <- 7
  g.p$sdd.rate <- 1.4
}

# load landscape
lc.df <- read.csv(paste0("data/gb/spread_", res, ".csv"))
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)

# initialize
cell.init <- get_pt_id(lc.df, c(739235.9, 4753487)) # 1922: herbarium_records.R


########
## Run model
########
# run sensitivity analysis
out.dir <- paste0("out/sdd_tune/", res, "/4_5/")
global_sensitivity(par.ls, nSamp, ngrid, ncell, g.p, select(lc.df, -MinObsYear), 
                   sdd=NULL, cell.init, control.p=NULL, save_yrs=1:g.p$tmax,
                   verbose=T, sim.dir=paste0(out.dir, "sims/"))


########
## validate pattern of spread 
########

f <- dir(out.dir, "N", recursive=TRUE)
Years <- sort(unique(lc.df$MinObsYear)[!is.na(unique(lc.df$MinObsYear))])
k.hist.occ <- map(Years, ~which(lc.df$MinObsYear <= .))
if(!dir.exists(paste0(out.dir, "pred"))) dir.create(paste0(out.dir, "pred"))

p.c <- makeCluster(g.p$n.cores); registerDoSNOW(p.c)
foreach(i=seq_along(f), .packages=Packages) %dopar% {
  N <- readRDS(paste0(out.dir, f[i])) 
  par <- read.csv(paste0(out.dir, str_replace(str_replace(f[i], "N", "results"),
                                              ".rds", ".csv")))
  results <- merge(data.frame(Year=Years), par) %>%
    mutate(pCorr=map2_dbl(k.hist.occ, Years, ~sum(N[.x,.y-1922]>0)/length(.x)))
  write_csv(results, paste0(out.dir, "pred/results_", i, ".csv"))
  return(i)
}
stopCluster(p.c)

out <- map_dfr(dir(paste0(out.dir, "pred"), full.names=T), read.csv)
write_csv(out, paste0(out.dir, "dispersal_out.csv"))

if(plots) {
  library(viridis)
  # univariate
  ggplot(out, aes(sdd.rate, pCorr)) + 
    geom_vline(xintercept=g.p$sdd.rate) + 
    geom_point(alpha=0.1) + stat_smooth(method="loess", se=F) + 
    facet_wrap(~Year) + ylim(0,1) +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.rate and prediction of historical records")
  ggplot(out, aes(sdd.max, pCorr)) + 
    geom_boxplot(aes(group=sdd.max), outlier.size=0.5) +
    geom_vline(xintercept=g.p$sdd.max) + 
    stat_smooth(method="loess") + facet_wrap(~Year) + ylim(0,1) +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.max and prediction of historical records")
  ggplot(out, aes(n.ldd, pCorr)) + 
    geom_boxplot(aes(group=n.ldd), outlier.size=0.5) +
    geom_vline(xintercept=g.p$n.ldd) + 
    stat_smooth(method="loess") + facet_wrap(~Year) + ylim(0,1) +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("n.ldd and prediction of historical records")
  
  # bivariate
  ggplot(out, aes(sdd.rate, pCorr, colour=n.ldd, group=n.ldd)) + 
    geom_vline(xintercept=g.p$sdd.rate) + ylim(0,1) + 
    geom_point(alpha=0.4) + stat_smooth(method="loess", se=F) + 
    facet_wrap(~Year) + scale_colour_viridis(option="B") +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.rate, n.ldd, and prediction of historical records")
  ggplot(out, aes(sdd.max, pCorr, colour=n.ldd, group=n.ldd)) + 
    geom_vline(xintercept=g.p$sdd.max) + ylim(0,1) + 
    geom_point(alpha=0.4) + stat_smooth(method="loess", se=F) + 
    facet_wrap(~Year) + scale_colour_viridis(option="B") +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.max, n.ldd, and prediction of historical records")
  ggplot(out, aes(sdd.rate, pCorr, colour=sdd.max, group=sdd.max)) + 
    geom_vline(xintercept=g.p$sdd.rate) + ylim(0,1) + 
    geom_point(alpha=0.4) + stat_smooth(method="loess", se=F) + 
    facet_wrap(~Year) + scale_colour_viridis(option="B") +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.rate, sdd.max, and prediction of historical records")
  ggplot(out, aes(n.ldd, pCorr, colour=sdd.max, group=sdd.max)) + 
    geom_vline(xintercept=g.p$n.ldd) + ylim(0,1) + 
    geom_point(alpha=0.4) + stat_smooth(method="loess", se=F) + 
    facet_wrap(~Year) + scale_colour_viridis(option="B") +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("n.ldd, sdd.max, and prediction of historical records")
  
  # trivariate
  ggplot(out, aes(sdd.rate, pCorr, colour=sdd.max, group=sdd.max)) + 
    geom_vline(xintercept=g.p$sdd.rate) + ylim(0,1) + 
    geom_point(alpha=0.4) +  
    facet_grid(n.ldd~Year) + scale_colour_viridis(option="B") +
    labs(y="Proportion of correctly predicted presences") +
    ggtitle("sdd.rate, sdd.max, n.ldd, and prediction of historical records")
  
  # subset: final year
  ggplot(filter(out, Year==2018), aes(sdd.rate, pOcc, colour=n.ldd, group=n.ldd)) + 
    geom_point(alpha=0.8) + stat_smooth(method="loess", se=F) + 
    scale_colour_viridis(option="B")
  ggplot(filter(out, Year==2018), aes(sdd.max, pOcc, colour=n.ldd, group=n.ldd)) + 
    geom_point(alpha=0.8) + stat_smooth(method="loess", se=F) + 
    scale_colour_viridis(option="B")
  ggplot(filter(out, Year==2018), aes(n.ldd, pOcc)) + geom_boxplot(aes(group=n.ldd)) +
    stat_smooth(method="loess")
  
  
  ggplot(filter(out, Year==2018), aes(sdd.rate, pCorr, colour=n.ldd, group=n.ldd)) + 
    geom_point(alpha=0.8) + stat_smooth(method="loess", se=F, size=1.5) + 
    scale_colour_viridis(option="B")
  ggplot(filter(out, n.ldd==15), aes(sdd.rate, pCorr, colour=sdd.max, group=sdd.max)) + 
    geom_point(alpha=0.8) + 
    scale_colour_viridis(option="B")
  ggplot(filter(out, Year==2018), aes(sdd.rate, pCorr, colour=sdd.max)) + 
    geom_point(alpha=0.8) + facet_wrap(~n.ldd) + ylim(0,1) +
    scale_colour_viridis(option="B")
  ggplot(filter(out, Year==2018), aes(sdd.rate, pCorr, colour=n.ldd, group=n.ldd)) + 
    geom_point(alpha=0.8) + facet_wrap(~sdd.max) + ylim(0,1) +
    scale_colour_viridis(option="B")
  
  ggplot(out, aes(sdd.rate, pCorr, colour=Year, group=Year)) + 
    geom_point(alpha=0.8) +
    facet_wrap(~n.ldd) +
    stat_smooth(method="loess", se=F)
  
  out %>% filter(pCorr > 0.8 & Year==2018) %>% arrange(pCorr) %>% summary()
  # medians:
  # sdd.rate = 0.0589
  # sdd.max = 27
  # n.ldd = 13
  # pOcc = 0.671
  # pSB = 0.814
  # pK = 0.537
  # pCorr = 0.818
}



# N.df <- map(dir(paste0(out.dir, "sims"), "N", full.names=T), 
#             ~data.frame(N=c(readRDS(.)), 
#                         Year=rep(1:g.p$tmax, each=ngrid))) %>%
#   do.call("rbind", .) %>%
#   mutate(id=rep(1:ngrid, times=length(dir(paste0(out.dir, "sims"), "N"))*g.p$tmax),
#          Occ=N>0) %>%
#   group_by(id, Year) %>%
#   summarise(pOcc=mean(Occ), mnN=mean(N))


