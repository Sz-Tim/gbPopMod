# USDA-NIFA Buckthorn model
# Tim Szewczyk

# This is an example of how the buckthorn model could be called in the economic
# model. Part I assumes the model is written in R (or can directly call an R
# function) and iterates the buckthorn population one time step, alternating
# with landowner interactions and decisions. Part II assumes the model is
# written in another language and the R script is called, with inputs and
# outputs written to external files in each iteration.

# The buckthorn model functions are stored as an R package called gbPopMod
# hosted on GitHub. Prior to publication, the repository is private. You can
# install the package along with all other required packages with:
if(!requite(gbPopMod)){
  devtools::install_github("Sz-Tim/gbPopMod", dependencies=TRUE,
                         auth_token="886b37e1694782d91c33da014d201a55d0c80bfb")
}
# The following packages are called by various gbPopMod functions:
# - here: easier file directory navigation
# - doSNOW: running in parallel
# - foreach: running in parallel
# - fastmatch: faster SDD neighborhood calculation
# - tidyverse: tidyr, dplyr, ggplot, purrr, tibble, stringr, forcats, readr
# - magrittr: additional pipe (%>%) functions
# - *scales: generating color scales in some plotting functions (PLOTS ONLY)
# - *gganimate: generating gifs (PLOTS ONLY)
# - *dismo: boosted regression trees for sensitivity analysis (SENSITIVITY ONLY)
Packages <- c("here", "gbPopMod", "tidyverse", "magrittr")
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
lc.df <- suppressMessages(read_csv("data/20a_full.csv")) # test: 9km_car.csv; full: 20a_full.csv
ngrid <- nrow(lc.df)
ncell <- sum(lc.df$inbd)
id.i <- lc.df %>% select(id, id.in)
# Note that lc.df and id.i contain two id columns:
# - id: rectangular grid ID; identified for all cells
# - id.in: inbound ID; identified for inbound cells, NA if out of bounds
# This is to allow simpler spatial calculations (e.g., SDD neighborhoods)

#--- set parameters as default; ?set_control_p; ?set_g_p
tmax <- 3
g.p <- set_g_p(m=c(7,7,3,3,3,3), g.D=0)
c.p <- set_control_p(null_ctrl=FALSE, 
                     pTrt.man=0.2,  # 5% of cells cut and/or spray
                     pTrt.grd=0.2,  # 5% of cells use ground cover
                     lc.chg=TRUE,
                     pChg=.01  # 1% of cells harvest forest
)

#--- calculate SDD neighborhoods & initialize buckthorn
sdd <- sdd_set_probs(ncell, lc.df, g.p, verbose=T)
N.init <- pop_init(ngrid, g.p, lc.df)  # eventually will be stored & loaded
B <- matrix(0, nrow=ngrid, ncol=tmax+1) # [cell, year]
N <- array(0, dim=c(ngrid, tmax+1, 6, max(g.p$m))) # [cell, year, LC, age]
N[,1,,] <- N.init
for(l in 1:6) { if(g.p$m[l] < 7) { N[,,l,g.p$m[l]:6] <- NA } }


for(t in 1:tmax) {
  ##---
  ## 2. Management plans are decided
  ##---
  # economic decision model
  # - which cells perform which management actions?
  # - here, management is randomized for illustration
  # Instead of using trt_assign(), the economic model can generate a two-column
  # dataframe where id (= id.i$id) and Trt (= treatment method).
  # Forest harvest is slightly more complicated, requiring a dataframe with the
  # id and id.in for the cells harvesting, and a dataframe with the change in 
  # each forest type
  
  #--- harvest timber
  # identify which cells will harvest, and how much of each forst type
  f_cut.i <- cut_assign(pChg=c.p$pChg, ncell=ncell, lc.df=lc.df, forest.col=6:9)
  # update lc.df with new forest proportions
  lc.df[f_cut.i$id.chg$id,] <- cut_forest(f_cut.i$id.chg, f_cut.i$mx, 
                                          forest.col=6:9, lc.df)
  # update SDD probabilities based on bird preferences for new LC composition
  sdd.i <- tibble(id.in=unique(
    arrayInd(which(sdd$i %in% f_cut.i$id.chg$id.in), dim(sdd$i))[,4]), 
    id=id.i$id[match(id.in, id.i$id.in)])
  sdd_new <- sdd_set_probs(nrow(sdd.i), lc.df, g.p, sdd.i)
  sdd$i[,,,sdd.i$id.in] <- sdd_new$i
  sdd$sp[sdd.i$id.in] <- sdd_new$sp
  
  #--- identify cells for groundcover treatments
  grd_cover.i <- trt_assign(id.i=id.i, ncell=ncell, 
                            pTrt=c.p$pTrt.grd, trt.eff=c.p$grd.trt)
  
  #--- identify cells for cut/spray treatments
  cut_spray.i <- trt_assign(id.i=id.i, ncell=ncell, 
                            pTrt=c.p$pTrt.man, trt.eff=c.p$man.trt)
  
  
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
  out <- iterate_pop(ngrid, ncell, N[,t,,], B[,t], g.p, lc.df, sdd, control.p, 
                     grd_cover.i, cut_spray.i, read_write=FALSE, path=NULL)
  N[,t+1,,] <- out$N
  B[,t+1] <- out$B
}

# visualize output
N.tot <- apply(N[,,,7], 1:2, sum) # sum adults across land cover categories
out.df <- lc.df %>%
  mutate(N.0=N.tot[,1],
         B.0=B[,1],
         N.final=N.tot[,tmax+1],
         B.final=B[,tmax+1])
# cell abundances through time
matplot(t(N.tot), type="l", lty=1, col=rgb(0,0,0,0.3))
matplot(t(B), type="l", lty=1, col=rgb(0,0,0,0.3))
# final maps
theme_set(theme_bw())
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=N.final)) + geom_point(aes(colour=N.0>0)) +
  scale_fill_gradient(low="white", high="red") +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue"))
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=B.final)) + geom_point(aes(colour=N.0>0)) +
  scale_fill_gradient(low="white", high="red") +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue"))
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=N.final>0)) + geom_point(aes(colour=N.0>0)) +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue")) +
  scale_fill_manual(values=c("FALSE"="gray80", "TRUE"="red"))
ggplot(out.df, aes(x=lon, y=lat)) + 
  geom_tile(aes(fill=B.final>0)) + geom_point(aes(colour=N.0>0)) +
  scale_colour_manual(values=c("FALSE"=NA, "TRUE"="blue")) +
  scale_fill_manual(values=c("FALSE"="gray80", "TRUE"="red"))

