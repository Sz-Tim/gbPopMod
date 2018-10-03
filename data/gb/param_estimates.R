# Estimation of demographic parameters
# Tim Szewczyk

# For the USDA glossy buckthorn project, parameters are either global or vary
# with land cover type. These varying parameters are assigned with a vector, 
# where the elements correspond with the following land cover categories:
# - Open 
# - Other (non-suitable)
# - Deciduous Forest 
# - Other Evergreen Forest (excludes White Pine)
# - White Pine Forest
# - Mixed Forest
# The plausible range and best estimate for each parameter are informed by
# either field experiments, lab experiments, published literature on glossy
# buckthorn, published literature on similar species, expert opinion, or some
# combination.


########
## set up
########
library(tidyverse); library(magrittr); library(gbPopMod)
res <- c("20ac", "9km2")[2]
# Data sets
fec_allen <- read_csv("data/gb/Allen_fecundity.csv") %>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
flwr_allen <- read_csv("data/gb/Allen_flower.csv")
dens_lee <- read_csv("data/gb/Lee_density.csv") %>% group_by(Canopy)
fec_lee <- read_csv("data/gb/Lee_fecundity.csv")
emerge_lee <- read_csv("data/gb/Lee_emerge.csv")
germ_lee <- read_csv("data/gb/Lee_germ.csv")
ann_aiello <- read_csv("data/gb/Aiello_fral-df-ann.csv") %>%
  mutate(nhlc=factor(nhlc, labels=c("Opn", "Dec", "Mxd", "WP"))) %>%
  filter(!is.na(nhlc))
edd_maps <- read.csv("data/gb/eddmaps_gb.csv", 
                     na.strings=c("NA", "NULL", "")) %>%
  filter(grossAreaInAcres > infestedAreaInAcres &
           grossAreaInAcres > 15 & grossAreaInAcres < 500) %>%
  mutate(prop_infested=infestedAreaInAcres/grossAreaInAcres)
# Storage
par.best <- set_g_p()
par.rng <- set_sensitivity_pars(names(par.best))


########
## Fecundity parameters
########

## p.f: Mean flowering probability
#--- Source: Allen Lab field + Aiello-Lammens field
p.f_allen <- flwr_allen %>% group_by(Plot_type) %>% summarise(p.f=mean(Fruit))
p.f_aiello <- ann_aiello %>% group_by(nhlc) %>% summarise(p.f=mean(fruit>0))
par.best$p.f <- c(mean(c(p.f_allen$p.f[1], p.f_aiello$p.f[1])), 0,
                  p.f_aiello$p.f[2], p.f_allen$p.f[3], p.f_allen$p.f[3],
                  mean(c(p.f_allen$p.f[2], p.f_aiello$p.f[3])))
par.rng$p.f$min <- par.best$p.f * 0.75
par.rng$p.f$max <- pmin(1, par.best$p.f * 1.25)


## mu: Mean number of fruits from a flowering adult
#--- Source: Tom Lee field + Allen Lab field + Aiello-Lammens field
# issues: 
#  - Lee = cleared only (ideal conditions)
#  - Allen = only 2 branches; cleared + mixed + wp
#  - Aiello-Lammens = cleared, deciduous, mixed
#  - No data for Other, Evg
prop <- c("open"=NA, "closed"=NA)
mu_lee <- filter(fec_lee, nFruit>0)$nFruit
mu_aiello <- filter(ann_aiello, fruit>0) %>% group_by(nhlc) %>% 
  summarise(mu=mean(fruit))
mu_allen <- fec_allen %>% filter(Treatment=="closed") %>%
  group_by(Plot_type, Plant) %>%
  summarise(tot_fruit=sum(T_t1, na.rm=TRUE)) %>% ungroup %>%
  mutate(canopy=case_when(Plot_type == "cleared" ~ "open",
                          Plot_type != "cleared" ~ "closed"))
prop[1] <- mean(filter(mu_allen, Plot_type=="cleared")$tot_fruit)/mean(mu_lee)
prop[2] <- mean(filter(mu_allen, Plot_type=="mixed")$tot_fruit)/mu_aiello$mu[3]
mu_allen <- mu_allen %>% 
  mutate(est_fruit=case_when(canopy == "open" ~ tot_fruit/prop[1],
                             canopy == "closed" ~ tot_fruit/prop[2])) %>% 
  group_by(Plot_type) %>% summarise(mn=mean(est_fruit))
par.best$mu <- c(mean(mu_lee), 0, mu_aiello$mu[2], 
                 mu_allen$mn[3], mu_allen$mn[3], mu_allen$mn[2])
par.rng$mu$min <- par.best$mu * 0.75
par.rng$mu$max <- par.best$mu * 1.25
  

## gamma: Mean number of seeds per fruit
#--- Source: Tom Lee field
frt_data <- filter(fec_lee, nFruit>0)$nSeedFruit
par.best$gamma <- mean(frt_data)
par.rng$gamma$min <- quantile(frt_data, 0.25)
par.rng$gamma$max <- quantile(frt_data, 0.75)


## m: Age at adulthood
#--- Source: Expert opinion (Tom Lee)
par.best$m <- c(3, 3, 7, 7, 7, 7)
par.rng$m$min <- c(2, 2, 4, 4, 4, 4)
par.rng$m$max <- c(4, 4, 8, 8, 8, 8)




########
## Dispersal parameters
########

## p.c: Proportion of fruits consumed by birds
#--- Source: Allen Lab field -- see data/gb/pc_calculation.R
par.best$p.c <- c(0.149, 0.273, 0.233)[c(1,1,2,3,3,2)]
par.rng$p.c$min <- c(0.113, 0.222, 0.205)[c(1,1,2,3,3,2)]
par.rng$p.c$max <- c(0.206, 0.314, 0.25)[c(1,1,2,3,3,2)]


## sdd.rate: Rate parameter for SDD exponential kernel (unit = cells)
#--- Source: Merow et al 2011 (mean distance = 2.14km)
par.best$sdd.rate <- ifelse(res=="20ac", 
                            1/(2.14/0.285),  # 20 acre cell ~ 0.285 x 0.285 km
                            1/(2.14/3))  # 9 km^2 cell = 3 x 3 km
par.rng$sdd.rate$min <- 1/(1/par.best$sdd.rate * 1.25)
par.rng$sdd.rate$max <- 1/(1/par.best$sdd.rate * 0.75)


## sdd.max: Maximum SDD distance
#--- Source: Merow et al 2011 (max distance ~ 22km, but )
par.best$sdd.max <- ifelse(res=="20ac", 
                           floor(22/0.285),  # 20 acre cell ~ 0.285 x 0.285 km
                           floor(22/3))  # 9 km^2 cell = 3 x 3 km
par.rng$sdd.max$min <- floor(par.best$sdd.max * 0.75)
par.rng$sdd.max$max <- ceiling(par.best$sdd.max * 1.25)


## bird.hab: relative bird habitat preferences
#--- Source: Merow et al 2011 (Dev=.39, Ag=.44, Dec=.06, Evg=.11)
par.best$bird.hab <- c(.39, .44, .06, .11, .11, .11)
par.best$bird.hab <- par.best$bird.hab/sum(par.best$bird.hab)
par.rng$bird.hab$min <- (par.best$bird.hab*0.75)
par.rng$bird.hab$max <- (par.best$bird.hab*1.25)


## n.ldd: Annual number of long distance dispersal events
#--- Source: Merow et al 2011 (1)
par.best$n.ldd <- 1
par.rng$n.ldd$min <- 0
par.rng$n.ldd$max <- 5



########
## Survival parameters
########

## s.c: Survival of consumed seeds
#--- Source: Bartuszevige & Gorchov 2006, LaFleur et al 2009, Ramirez & Ornelas 2009
s.c_lit <- c(0.86, 0.47, 0.56, 0.54, 0.46, 0.62)
par.best$s.c <- mean(s.c_lit)
par.rng$s.c$min <- quantile(s.c_lit, 0.25)
par.rng$s.c$max <- quantile(s.c_lit, 0.75)


## s.B: Annual survival rate in seed bank
#--- Source: Expert opinion + CITATION
# life expectancy = 3 years
par.best$s.B <- exp(-1/3)  # life expectancy = -1/log(annual survival)
par.rng$s.B$min <- exp(-1/(3*0.75))
par.rng$s.B$max <- exp(-1/(3*1.25))


## s.M: Annual survival rate for juveniles
#--- Source: Literature (Ranges from similar species)
par.rng$s.M$min <- par.best$s.M * 0.75
par.rng$s.M$max <- pmin(par.best$s.M * 1.25, 1)


## s.N: Annual survival rate for adults
#--- Source: Expert opinion -- no observed mortality in USDA time frame
par.best$s.N <- rep(1, 6)
par.rng$s.N$min <- rep(0.9, 6)
par.rng$s.N$max <- rep(1, 6)


## K: Carrying capacity for adults
#--- Source: Tom Lee field, EDDMapS
em <- max(edd_maps$prop_infested)
dens_data <- dens_lee %>% summarise(mn=mean(n_g1m_ha),
                                    q25=quantile(n_g1m_ha, 0.25),
                                    q75=quantile(n_g1m_ha, 0.75))
dens_data[,2:4] <- dens_data[,2:4] * ifelse(res=="20ac", 8.1, 900) * em
par.best$K <- with(dens_data, c(mn[2], 0, mn[1], mn[1], mn[1], mn[1]))
par.rng$K$min <- with(dens_data, c(q25[2], 0, q25[1], q25[1], q25[1], q25[1]))
par.rng$K$max <- with(dens_data, c(q75[2], 100, q75[1], q75[1], q75[1], q75[1]))



########
## Seedling parameters
########

## g.D: Proportion of seeds germinating in same year they were produced
#--- Source:
par.best$g.D <- 0
par.rng$g.D$min <- 0
par.rng$g.D$max <- 0


## g.B: Proportion of seeds germinating from the seed bank
#--- Source: Tom Lee lab
germ_data <- filter(germ_lee, Treatment!="ColdTemp")$propGerm
par.best$g.B <- mean(germ_data)
par.rng$g.B$min <- quantile(germ_data, 0.25)
par.rng$g.B$max <- quantile(germ_data, 0.75)


## p: Proportion of germinants that establish
#--- Source: Tom Lee field + Allen Lab field
par.best$p <- c(0.07, 0, 0.08, 0.02, 0.02, 0.03)/par.best$g.B
par.rng$p$min <- par.best$p * 0.75
par.rng$p$max <- pmin(1, par.best$p * 1.25)






