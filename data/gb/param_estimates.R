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
res <- c("20ac", "9km2")[1]
# Data sets
fec_allen <- read_csv("data/gb/Allen_fecundity.csv") %>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
flwr_allen <- read_csv("data/gb/Allen_flower.csv")
dens_lee <- read_csv("data/gb/Lee_density.csv") %>% group_by(Canopy)
fec_lee <- read_csv("data/gb/Lee_fecundity.csv")
emerge_lee <- read_csv("data/gb/Lee_emerge.csv")
germ_lee <- read_csv("data/gb/Lee_germ.csv")
sM_lee <- read_csv("data/gb/Lee_juvMortality.csv")
flwr_bibaud <- read_csv("data/gb/Bibaud_flower.csv")
ann_aiello <- read_csv("data/gb/Aiello_fral-df-ann.csv") %>%
  mutate(nhlc=factor(nhlc, labels=c("Opn", "Dec", "Mxd", "WP"))) %>%
  filter(!is.na(nhlc))
edd_maps <- read.csv("data/gb/eddmaps_gb.csv", 
                     na.strings=c("NA", "NULL", "")) %>%
  filter(grossAreaInAcres > infestedAreaInAcres &
           grossAreaInAcres > 15 & grossAreaInAcres < 500) %>%
  mutate(prop_infested=infestedAreaInAcres/grossAreaInAcres)
mort_eis <- read_csv("data/gb/Eisenhaure_cut.csv") %>% 
  mutate(Date=as.Date(Date, format="%m/%d/%Y"),
         Season=forcats::lvls_reorder(Season, c(3,2,1)))
# Storage
par.best <- set_g_p()
par.rng <- set_sensitivity_pars(names(par.best))


########
## Fecundity parameters
########

## p.f: Mean flowering probability
#--- Source: Allen Lab field + Aiello-Lammens field + Bibaud field
p.f_allen <- flwr_allen %>% group_by(Plot_type) %>% summarise(p.f=mean(Fruit))
p.f_aiello <- ann_aiello %>% group_by(nhlc) %>% summarise(p.f=mean(fruit>0))
m_bibaud <- flwr_bibaud %>% filter(!is.na(Flowering)) %>% summarise(m=min(Age))
p.f_bibaud <- flwr_bibaud %>% filter(Age > m_bibaud$m) %>% group_by(Canopy) %>% 
  summarise(nFlower=sum(Flowering, na.rm=T),
            nNotFlower=sum(NotFlowering, na.rm=T),
            p.f=nFlower/nNotFlower)
par.best$p.f <- c(mean(c(p.f_allen$p.f[1], p.f_aiello$p.f[1])), 
                  0,
                  mean(p.f_aiello$p.f[2], p.f_bibaud$p.f[1]), 
                  p.f_bibaud$p.f[2], p.f_bibaud$p.f[2],
                  mean(c(p.f_bibaud$p.f, p.f_aiello$p.f[3])))
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
#--- Source: Bibaud field + Expert opinion (Tom Lee)
m_bibaud <- flwr_bibaud %>% filter(!is.na(Flowering)) %>% summarise(m=min(Age))
par.best$m <- c(3, 3, 7, 7, 7, 7)
par.rng$m$min <- c(2, 2, 4, 4, 4, 4)
par.rng$m$max <- c(4, 4, 10, 10, 10, 10)




########
## Dispersal parameters
########

## p.c: Proportion of fruits consumed by birds
#--- Source: Allen Lab field -- see data/gb/pc_calculation.R
par.best$p.c <- c(0.165, 0.296, 0.252)[c(1,1,2,3,3,2)]
par.rng$p.c$min <- c(0.122, 0.25, 0.219)[c(1,1,2,3,3,2)]
par.rng$p.c$max <- c(0.229, 0.333, 0.262)[c(1,1,2,3,3,2)]


## sdd.rate: Rate parameter for SDD exponential kernel (unit = cells)
#--- Source: Merow et al 2011 (mean distance = 2.14km)
par.best$sdd.rate <- ifelse(res=="20ac", 
                            1/(2.14/0.285),  # 20 acre cell ~ 0.285 x 0.285 km
                            1/(2.14/3))  # 9 km^2 cell = 3 x 3 km
par.rng$sdd.rate$min <- 1/(1/par.best$sdd.rate * 1.25)
par.rng$sdd.rate$max <- 1/(1/par.best$sdd.rate * 0.75)


## sdd.max: Maximum SDD distance
#--- Source: Merow et al 2011 (max distance ~ 22km, but seems unnecessary here)
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
p_lee <- emerge_lee %>% group_by(Treatment) %>% summarise(p=mean(propEmerge))
par.best$p <- c(p_lee$p[2], 0, p_lee$p[3], 
                p_lee$p[4], p_lee$p[4], mean(p_lee$p[3:4])) / par.best$g.B
par.rng$p$min <- par.best$p * 0.75
par.rng$p$max <- pmin(1, par.best$p * 1.25)
# ground cover treatment: cover crop
p_lee$p[5]/par.best$g.B



########
## Management mortality
########

mort_none <- mort_eis %>% 
  filter(Season == "Start" & Treatment=="NONE") %>%
  mutate(Y0=case_when(Year=="2008" ~ 2008,
                      Year=="2009" ~ 2009,
                      Year=="2010" ~ 2009)) %>%
  bind_rows(., filter(., Year==2009) %>% mutate(Y0=2008)) %>%
  mutate(category=case_when(Y0==Year & Height=="short" ~ "S0",
                            Y0<Year & Height=="short" ~ "S1",
                            Y0==Year & Height=="tall" ~ "T0",
                            Y0<Year & Height=="tall" ~ "T1")) %>%
  select(Y0, Subblock, Treatment, N, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0)
mort_cut <- mort_eis %>% 
  filter(Season == "Start" & Treatment=="CUT") %>%
  mutate(Y0=case_when(Year=="2008" ~ 2008,
                      Year=="2009" ~ 2009,
                      Year=="2010" ~ 2009)) %>%
  bind_rows(., filter(., Year==2009) %>% mutate(Y0=2008)) %>%
  mutate(category=case_when(Y0==Year & Height=="short" ~ "S0",
                            Y0<Year & Height=="short" ~ "S1",
                            Y0==Year & Height=="tall" ~ "T0",
                            Y0<Year & Height=="tall" ~ "T1")) %>%
  select(Y0, Subblock, Treatment, N, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0)
mean((mort_cut$pN - mort_none$pN)) # 97% decrease overall, but one big outlier
mean((mort_cut$pN - mort_none$pN)[1:3]) # 68% decrease for Y0=2008
mean((mort_cut$pN - mort_none$pN)[4:6]) # 125% decrease for Y0=2009



stop("Unused exploration beyond this point.")
##------------- cut efficacy exploration -------------##
## Notation:
##  S_: 'short' abundance
##  T_: 'tall' abundance
##  N_: total abundance
##  _0: start of year 0
##  _1: start of year 1
##  s: survival rate (no treatment)
##  g: pr(S0 -> T1 | survival)
##  B: new sprout abundance
##  r: pr(resprout from cutting) = 1 - mortality from cutting
##  p: proportion of resprouts from T0 to S1

## No Treatment
# N0 = S0 + T0
# N1 = S1 + T1 = N0*s + B
# S1 = S0*s - S0*s*g + B
# T1 = T0*s + S0*s*g
#   g = (T1-T0*s)/(S0*s)
#   B = N1 - N0*s

## Cut
# N0 = S0 + T0
# N1 = S1 + T1 = N0*s*r + B
# S1 = S0*s*r - S0*s*r*g + T0*s*r*p + B
# T1 = T0*s*r*(1-p) + S0*s*r*g
#   r = (N1-B)/(N0*s)
#   p = 1 - (T1-S0*s*r*g)/(T0*s*r)

## M: mechanical only
# Estimate g, B from Treatment=="NONE" for each subblock
s <- 0.6  # need s ≤ 0.18 for all g>0 & B>0
mort_height <- mort_eis %>% 
  filter(Season == "Start" & Treatment=="NONE") %>%
  mutate(Y0=case_when(Year=="2008" ~ 2008,
                      Year=="2009" ~ 2009,
                      Year=="2010" ~ 2009)) %>%
  bind_rows(., filter(., Year==2009) %>% mutate(Y0=2008)) %>%
  mutate(category=case_when(Y0==Year & Height=="short" ~ "S0",
                            Y0<Year & Height=="short" ~ "S1",
                            Y0==Year & Height=="tall" ~ "T0",
                            Y0<Year & Height=="tall" ~ "T1")) %>%
  select(Y0, Subblock, Treatment, N, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0,
         g=(T1-T0*s)/(S0*s),
         B=N1-N0*s, pB=B/N0)
none_est <- mort_height %>% 
  group_by(Treatment) %>% 
  summarise(mn_g=mean(g), mn_pB=mean(pB))

# Calculate r, p using estimates of g, B
mort_cut <- mort_eis %>% 
  filter(Season == "Start" & Treatment=="CUT") %>%
  mutate(Y0=case_when(Year=="2008" ~ 2008,
                      Year=="2009" ~ 2009,
                      Year=="2010" ~ 2009)) %>%
  bind_rows(., filter(., Year==2009) %>% mutate(Y0=2008)) %>%
  mutate(category=case_when(Y0==Year & Height=="short" ~ "S0",
                            Y0<Year & Height=="short" ~ "S1",
                            Y0==Year & Height=="tall" ~ "T0",
                            Y0<Year & Height=="tall" ~ "T1")) %>%
  select(Y0, Subblock, Treatment, N, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0,
         g=none_est$mn_g[match(Subblock, none_est$Subblock)],
         pB=none_est$mn_pB[match(Subblock, none_est$Subblock)],
         B=pB*N0,
         r=(N1-B)/(N0*s),
         p=1-(T1-S0*s*r*g)/(T0*s*r))
mean((mort_cut$pN - mort_height$pN)[1:3])

# WITHIN YEAR
mort_height <- mort_eis %>% 
  filter(Season != "Mid" & Year !=2010 & Treatment=="NONE") %>%
  mutate(category=case_when(Season=="Start" & Height=="short" ~ "S0",
                            Season=="End" & Height=="short" ~ "S1",
                            Season=="Start" & Height=="tall" ~ "T0",
                            Season=="End" & Height=="tall" ~ "T1")) %>%
  select(Subblock, Treatment, N, Year, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0,
         g=(T1-T0*s)/(S0*s),
         B=N1-N0*s, pB=B/N0)
none_est <- mort_height %>% 
  group_by(Subblock, Treatment) %>% 
  summarise(mn_g=mean(g), mn_pB=mean(pB))

# Calculate r, p using estimates of g, B
mort_cut <- mort_eis %>% 
  filter(Season != "Mid" & Year !=2010 & Treatment=="CUT") %>%
  mutate(category=case_when(Season=="Start" & Height=="short" ~ "S0",
                            Season=="End" & Height=="short" ~ "S1",
                            Season=="Start" & Height=="tall" ~ "T0",
                            Season=="End" & Height=="tall" ~ "T1")) %>%
  select(Subblock, Treatment, N, Year, category) %>%
  spread(category, N) %>%
  mutate(N0=S0+T0, N1=S1+T1,
         dS=S1-S0, dT=T1-T0, dN=N1-N0,
         pS=dS/S0, pT=dT/T0, pN=dN/N0,
         g=none_est$mn_g[match(Subblock, none_est$Subblock)],
         pB=none_est$mn_pB[match(Subblock, none_est$Subblock)],
         B=pB*N0,
         r=(N1-B)/(s*N0),
         p=1-(T1-s*r*g*S0)/(s*r*T0))
mean(mort_cut$pN - mort_height$pN)

mort_sum <- mort_eis %>% group_by(Subblock, Treatment, Date, Season, Year) %>%
  summarise(N=sum(N))
mort_pct <- mort_eis %>% 
  filter(Season != "Mid" & Year != 2010) %>% 
  select(-Date) %>% 
  spread(Season, N) %>%
  mutate(chg.N=End-Start,
         chg.pct=chg.N/(Start))
mort_sum_pct <- mort_sum %>% ungroup %>%
  filter(Season != "Mid" & Year != 2010) %>% 
  select(-Date) %>% 
  spread(Season, N) %>%
  mutate(chg.N=End-Start,
         chg.pct=chg.N/(Start))
mort_diff <- mort_pct %>%
  select(Subblock, Treatment, Year, Height, chg.pct) %>%
  spread(Treatment, chg.pct) %>%
  mutate(chg.pct.cut_none=CUT-NONE)
mort_sum_diff <- mort_sum_pct %>%
  select(Subblock, Treatment, Year, chg.pct) %>%
  spread(Treatment, chg.pct) %>%
  mutate(chg.pct.cut_none=CUT-NONE)
