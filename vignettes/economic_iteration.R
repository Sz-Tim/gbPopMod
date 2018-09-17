# USDA-NIFA Buckthorn model
# Tim Szewczyk

# This is an example of how the buckthorn model could be called in the economic
# model. Part I assumes the model is written in R (or can directly call an R
# function) and iterates the buckthorn population one time step, alternating
# with landowner interactions and decisions. 

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
# devtools::install_github("Sz-Tim/gbPopMod", dependencies=T,
#                        auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
# help(package="gbPopMod")
# The following packages are called by various gbPopMod functions:
# - here: easier file directory navigation
# - doSNOW: running in parallel
# - foreach: running in parallel
# - fastmatch: faster SDD neighborhood calculation
# - tidyverse: tidyr, dplyr, ggplot, purrr, tibble, stringr, forcats, readr
# - magrittr: additional pipe (%>%) functions
# - *scales: generating color scales in some plotting functions (PLOTS ONLY)
# - *gganimate: generating gifs (PLOTS ONLY)
# - *viridis: for pretty colors (PLOTS ONLY)
# - *dismo: boosted regression trees for sensitivity analysis (SENSITIVITY ONLY)
Packages <- c("here", "gbPopMod", "tidyverse", "magrittr", "viridis")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))


########
## Part I: Iterate buckthorn model within R
########
# Here, the economic model is written in R and the buckthorn population is
# iterated each time step. The basic flow is:
# 1. Set up and initialization. Store initial buckthorn population, etc
# 2. Management plans are decided
# 3. Management is implemented, buckthorn grows and spreads
# 4. Repeat 2-3 for tmax years

##---
## 1. SET UP AND INITIALIZATION
##---
#--- load libraries & landscape
gifs <- TRUE
res <- "9km2"  # 9km2 or 20ac
load(paste0("data/USDA_", res, ".rda"))  # loads dataframe as lc.df
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)
id.i <- lc.df %>% select(id, id.in)
# Note that lc.df and id.i contain two id columns:
# - id: rectangular grid ID; identified for all cells
# - id.in: inbound ID; identified for inbound cells, NA if out of bounds
# This is to allow simpler spatial calculations (i.e., for SDD neighborhoods)

#--- set parameters
tmax <- 100
# cells that implement manual treatments
manual.i <- filter(lc.df, x>38 & x<45 & y>33 & y<40) %>% mutate(trt="m")
# cells that implement ground cover treatments
ground.i <- filter(lc.df, x>18 & x<25 & y>38 & y<45) %>% mutate(trt="g")
# cells that implement both
both.i <- filter(lc.df, x>28 & x<35 & y>28 & y<35) %>% mutate(trt="b")
mgmt.df <- rbind(manual.i, ground.i, both.i) %>% group_by(trt) %>%
  summarise(lon.mx=max(lon), lon.mn=min(lon),
            lat.mx=max(lat), lat.mn=min(lat))
mgmt.df <- with(mgmt.df, data.frame(trt=rep(trt, 4), 
                      lon=c(lon.mx, lon.mx, lon.mn, lon.mn),
                      lat=c(lat.mx, lat.mn, lat.mn, lat.mx)))
# global parameters: ?set_g_p
g.p <- set_g_p(tmax=tmax, n.cores=1,  
               K=c(5223180, 0, 770790, 770790, 770790, 770790),
               sdd.max=7, sdd.rate=1.4)
# control parameters: ?set_control_p
c.p <- set_control_p(null_ctrl=FALSE, 
                     t.trt=100,
                     grd.i=c(ground.i$id.in, both.i$id.in),
                     man.i=c(manual.i$id.in, both.i$id.in),
                     lc.chg=FALSE
)

#--- calculate SDD neighborhoods & initialize buckthorn
sdd <- sdd_set_probs(ncell, lc.df, g.p, verbose=T)
set.seed(1)
N.init <- pop_init(ngrid, g.p, lc.df, p.0=2000)
B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
N <- array(0, dim=c(ngrid, tmax+1, 6, max(g.p$m))) # [cell, year, LC, age]
N[,1,,] <- N.init
for(l in 1:6) { if(g.p$m[l] < 7) { N[,,l,g.p$m[l]:(max(g.p$m)-1)] <- NA } }
grd_cover.i <- mech_chem.i <- NULL


system.time({
for(t in 1:g.p$tmax) {
  ##---
  ## 2. Management plans are decided
  ##---
  # economic decision model
  # - which cells perform which management actions?
  # - here, management is randomized for illustration
  # Instead of using trt_assign(), the economic model can generate a two-column
  # dataframe where id (=id.i$id) and Trt (=treatment method).
  # Forest harvest is slightly more complicated, requiring a dataframe with the
  # id and id.in for the harvesting cells, and a dataframe with the change in 
  # each forest type
  
  #--- harvest timber
  if(c.p$lc.chg) {
    # identify which cells will harvest, and how much of each forst type
    f_cut.i <- cut_assign(pChg=c.p$pChg, ncell=ncell, 
                          lc.df=lc.df, forest.col=6:9)
    # update lc.df with new forest proportions
    lc.df[f_cut.i$id.chg$id,] <- cut_forest(f_cut.i$id.chg, f_cut.i$mx, 
                                            forest.col=6:9, lc.df)
    # update SDD probabilities based on bird preferences for new LC composition
    sdd.alter <- tibble(id.in=unique(
      arrayInd(which(sdd$i %in% f_cut.i$id.chg$id.in), dim(sdd$i))[,4]), 
      id=id.i$id[match(id.in, id.i$id.in)])
    sdd_new <- sdd_update_probs(lc.df, g.p, sdd.alter, sdd$i)
    sdd$i[,,1,sdd.alter$id.in] <- sdd_new$i
    sdd$sp[sdd.alter$id.in] <- sdd_new$sp
  }
  
  if(t >= c.p$t.trt) {
  #--- identify cells for groundcover treatments
  grd_cover.i <- trt_assign(id.i, assign_i=c.p$grd.i, trt.eff=c.p$grd.trt)
  
  #--- identify cells for cut/spray treatments
  mech_chem.i <- trt_assign(id.i, assign_i=c.p$man.i, trt.eff=c.p$man.trt)
  }
  
  
  ##---
  ## 3. MANAGEMENT IS IMPLEMENTED, BUCKTHORN GROWS & SPREADS
  ##---
  # buckthorn population model
  # Parameters are updated for cells that performed management actions, then
  # buckthorn grows and spreads. This is implemented as:
  # 1. Update parameters based on management actions
  # 2. Adults produce fruit
  # 3. Fruit is consumed and dispersed, or dropped
  # 4. Seeds germinate and establish
  g.p$n.ldd <- ifelse(t %% 5 == 0, 1, 0)
  out <- iterate_pop(ngrid, ncell, N[,t,,], B[,t], g.p, lc.df, sdd, c.p, 
                     grd_cover.i, mech_chem.i, read_write=FALSE, path=NULL)
  N[,t+1,,] <- out$N
  B[,t+1] <- out$B
  if(t %% 10 == 0) cat("Finished time", t, "\n")
}
})

# visualize output
plot.dir <- paste0("out/", res, "/econ/t_", g.p$tmax, "/trt_", c.p$t.trt)
if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive=T)
N.tot <- apply(N[,,,max(g.p$m)], 1:2, sum) # sum adults across LC categories
N.df <- as.data.frame(N.tot); names(N.df) <- 1:ncol(N.df)
out.all <- cbind(lc.df, N.df) %>%
  gather(year, N, 15:ncol(.)) %>% mutate(year=as.numeric(year), B=c(B))
out.df <- lc.df %>%
  mutate(N.0=N.tot[,1],
         B.0=B[,1],
         N.final=N.tot[,tmax+1],
         B.final=B[,tmax+1])
out.all <- filter(out.all, !is.na(lon))
write_csv(out.all, here(plot.dir, "out_all.csv"))
# cell abundances through time
theme_set(theme_bw() + theme(panel.background=element_rect(fill="gray30"),
                             panel.grid=element_blank()))
if(gifs) {
  library(gganimate)
  anim_save(here(plot.dir, "N.gif"),
            animate(ggplot(out.all, aes(x=lon, y=lat)) + 
                      geom_tile(aes(fill=log(N))) +
                      geom_polygon(data=mgmt.df, aes(colour=trt), fill=NA, size=1) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(paste("Treatment from\nyear", c.p$t.trt), 
                                          values=c("black", "blue", "purple"), 
                                          labels=c("both", "ground", "manual")) +
                      transition_time(year) +  
                      ggtitle("Adult abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
  anim_save(here(plot.dir, "B.gif"),
            animate(ggplot(out.all, aes(x=lon, y=lat)) + 
                      geom_tile(aes(fill=log(B))) +
                      geom_polygon(data=mgmt.df, aes(colour=trt), fill=NA, size=1) +
                      scale_fill_viridis(option="B") + 
                      scale_colour_manual(paste("Treatment from\nyear", c.p$t.trt), 
                                          values=c("black", "blue", "purple"), 
                                          labels=c("both", "ground", "manual")) +
                      transition_time(year) +  
                      ggtitle("Seed abundance. Year {frame_time}"),
                    nframes=n_distinct(out.all$year), 
                    width=800, height=600, units="px"))
}
# final maps
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=N.final)) + geom_point(aes(colour=N.0>0)) +
  scale_fill_viridis(option="B") +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue"))
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=log(B.final))) + geom_point(aes(colour=N.0>0)) +
  scale_fill_viridis(option="B") +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue"))
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=N.final>0)) + geom_point(aes(colour=N.0>0)) +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue")) +
  scale_fill_viridis(option="B", discrete=TRUE) 
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=B.final>0)) + geom_point(aes(colour=N.0>0)) +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue")) +
  scale_fill_viridis(option="B", discrete=TRUE)

