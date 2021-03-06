library(tidyverse)
fec_allen <- read.csv("data/gb/Allen_fecundity.csv") %>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))

# i_t1 = immature fruit at t=1
# m_t1 = mature fruit at t=1
# T_t1 = total fruit at t=1 (i+m)
# b_t1 = fruit in bag at t=1 (should accumulate)

# new fruit on closed branches:
#  on_branch_t1 = on_branch_t0 - new_in_bag + new_on_branch (+ consumed)
#  T_t1 = T_t0 - (b_t1-b_t0) + n_t1
#  n_t1 = T_t1 - T_t0 + b_t1 - b_t0
# ***Using total (vs. mature) fruit: 
#    - CW_WP Plant 4 has m=0 for all t
#    - some immature fruits might mature & drop between observations
# ***Sometimes n_t1 is negative -- assume underestimate of new_in_bag?

# new fruit on open branches:
#  on_branch_t1 = on_branch_t0 - fruit_dropped + new_on_branch - fruit_consumed
#  fruit_consumed = on_branch_t0 - on_branch_t1 + new_on_branch - fruit_dropped
#  c_e = T_t0 - T_t1 + n_e + d_e

# assume fruit production rate is constant among branches on each plant
#  new_rate = new_on_branch / (on_branch_t1 + new_in_bag + consumed)
#  n_rate = n_t1 / (T_t1 + (b_t1-b_t0) + c_e)
#  n_rate = n_e / (T_t1 + d_e + c_e)

# assume fruit drop rate is constant among branches on each plant
#  drop_rate = new_in_bag / (on_branch_t1 + new_in_bag + consumed)
#  d_rate = (b_t1-b_t0) / (T_t1 + (b_t1-b_t0) + c_e)
#  d_rate = d_e / (T_t1 + d_e + c_e)

# use plant-specific *_rate to estimate fruit_dropped, new_on_branch for open
#  n_e = n_rate*T_t1 + n_rate*d_e + n_rate*c_e
#  d_e = (d_rate*(T_t1+c_e))/(1-d_rate)
#  c_e = T_t0 - T_t1 + n_e + d_e
#      = T_t0 - T_t1 + 
#        n_rate*T_t1 + n_rate*(d_rate*(T_t1+c_e))/(1-d_rate) + n_rate*c_e + 
#        (d_rate*(T_t1+c_e))/(1-d_rate)
#      = (T_t0*(1-d_rate))/(1-n_rate) - T_t1


for(i in 2:n_distinct(fec_allen$t_1)) {
  for(j in 1:n_distinct(fec_allen$Plot_abbrev)) {
    for(k in 1:n_distinct(fec_allen$Plant)) {
      for(l in 1:n_distinct(fec_allen$Branch)) {
        i_t1 <- with(fec_allen, which(t_1==i & 
                                        Plot_abbrev==unique(Plot_abbrev)[j] &
                                        Plant==unique(Plant)[k] &
                                        Branch==unique(Branch)[l]))
        i_t0 <- with(fec_allen, which(t_1==(i-1) & 
                                        Plot_abbrev==unique(Plot_abbrev)[j] &
                                        Plant==unique(Plant)[k] &
                                        Branch==unique(Branch)[l]))
        fec_allen$i_t0[i_t1] <- fec_allen$i_t1[i_t0]
        fec_allen$m_t0[i_t1] <- fec_allen$m_t1[i_t0]
        fec_allen$T_t0[i_t1] <- fec_allen$T_t1[i_t0]
        fec_allen$b_t0[i_t1] <- fec_allen$b_t1[i_t0]
      }
    }
  }
}
fec_allen <- fec_allen %>%
  mutate(delta_T=T_t1-T_t0,
         delta_m=m_t1-T_t0,
         delta_b=b_t1-b_t0)
rates.df <- fec_allen %>%
  filter(Treatment=="closed") %>%
  mutate(n_t1=T_t1-T_t0+b_t1-b_t0) 
rates.df <- rates.df %>%
  mutate(delta_b=b_t1-b_t0) %>%
  group_by(Plot_abbrev, Plant) %>%
  summarise(n_rate=sum(n_t1, na.rm=T)/(sum(T_t1, delta_b, na.rm=T)),
            d_rate=sum(delta_b, na.rm=T)/(sum(T_t1, delta_b, na.rm=T))) %>%
  as.tibble

fec_all <- left_join(fec_allen, rates.df, by=c("Plot_abbrev", "Plant")) %>%
  filter(Treatment=="open") %>%
  mutate(c_e=(T_t0*(1-d_rate)/(1-n_rate) - T_t1),
         d_e=(d_rate*(T_t1+c_e))/(1-d_rate),
         n_e=n_rate*(T_t1 + d_e + c_e)) %>%
  group_by(Plot_type, Plot_abbrev, Plant) %>%
  summarise(ce=sum(c_e, na.rm=T),
            de=sum(d_e, na.rm=T),
            ne=sum(n_e, na.rm=T),
            p.c=sum(c_e, na.rm=T)/sum(T_t1, d_e, c_e, na.rm=T))
fec_all$p.c[fec_all$p.c<0] <- 0
fec_all$p.c[fec_all$p.c>1] <- 1
fec_all %>% ungroup %>% group_by(Plot_type) %>% 
  summarise(c_mn=mean(p.c, na.rm=T),
            c_md=median(p.c, na.rm=T),
            c_25=quantile(p.c, probs=0.25, na.rm=T),
            c_75=quantile(p.c, probs=0.75, na.rm=T))

























