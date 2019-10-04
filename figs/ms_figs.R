# Figures for the manuscript and talks
# Glossy buckthorn demographic CA model
# Tim Szewczyk

Packages <- c("tidyverse", "viridis", "gridExtra", "scales", "sf", "ggspatial")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
theme_set(theme_bw())

ms_fig.dir <- "../ms_DemogCA/figs/"
if(!dir.exists(ms_fig.dir)) dir.create(ms_fig.dir, recursive=T)
res <- c("20ac", "9km2")[1]
load(paste0("data/USDA_", res, ".rda"))

ms_fonts <- theme(plot.title=element_text(size=11),
                  axis.title=element_text(size=10),
                  axis.text=element_text(size=9),
                  legend.title=element_text(size=9),
                  legend.text=element_text(size=9))

# Supplements:
# - vignettes/Appendix_1.Rmd: Data and parameter details
# - vignettes/Appendix_2.Rmd: GSA and dispersal tuning
# - vignettes/Appendix_3.Rmd: Code to run models

# Figures:
# - mod_outline: model structure (produced in Keynote)
# - UNH_map: map of UNH properties, with New England inset
# - mod_output: maps of abundances, seed production, immigration, seed pressure
# - GSA_relInf: relative influence of parameters from GSA
# - mgmt_N_B: buckthorn N_adult & seed bank under 4 mgmt strategies
# - mgmt_propDiff: proportion difference compared to 'none'



########
## Appendixes
########
rmarkdown::render("vignettes/Appendix_1.Rmd", 
                  output_file="../../ms_DemogCA/supp/Appendix_1.pdf",
                  knit_root_dir="../")
rmarkdown::render("vignettes/Appendix_2.Rmd", 
                  output_file="../../ms_DemogCA/supp/Appendix_2.pdf",
                  knit_root_dir="../")
rmarkdown::render("vignettes/Appendix_3.Rmd", 
                  output_file="../../ms_DemogCA/supp/Appendix_3.pdf",
                  knit_root_dir="../")




########
## UNH_map
########
states <- maps::map("state", plot=F, fill=T) %>% st_as_sf(.) %>%
  filter(ID %in% c("maine", "new hampshire", "vermont", "massachusetts")) %>%
  st_set_crs(4326) %>% st_transform(crs=32618)
lc_unh.df <- read.csv(paste0("data/USDA_", res, "_mgmt.csv")) %>%
  mutate(N_adult=readRDS(paste0("data/inits/N_2018_", res, ".rds"))[.$id.full,6,7])
lc_unh.sf <- filter(lc_unh.df, !is.na(Property)) %>%
  st_as_sf(coords=c("lon", "lat")) %>%
  st_set_crs(32618) 
lc_bound <- dplyr::filter(lc.df, !is.na(lon)) %>% 
  dplyr::select(lon, lat, inbd) %>%
  raster::rasterFromXYZ() %>%
  raster::rasterToPolygons(dissolve=T) %>%
  fortify %>% 
  dplyr::select(long, lat, group)
lc_unh_bound <- filter(lc_unh.df, !is.na(lon)) %>% 
  select(lon, lat, inbd) %>%
  raster::rasterFromXYZ() %>%
  raster::rasterToPolygons(dissolve=T) %>%
  fortify %>% 
  select(long, lat, group)
lc_unh_bound.sf <- st_as_sf(lc_unh_bound, coords=c("long", "lat")) %>% 
  st_set_crs(32618)
towns.nh <- sf::st_read("data/UNH_woodlands/Towns/pbp.shp") %>%
  st_transform(crs=32618) %>% select(geometry)
towns.me <- sf::st_read("data/UNH_woodlands/Towns/Maine_Boundaries_Town_Lines.shp") %>%
  st_transform(crs=32618) %>% select(geometry)
towns <- rbind(towns.nh, towns.me) %>%
  st_crop(., lc_unh_bound.sf) 

context_map <- ggplotGrob(
  ggplot() + geom_sf(data=states, size=0.2) + 
    geom_polygon(data=lc_bound, aes(long, lat, group=group), 
                 fill=NA, colour="black") +
    geom_polygon(data=lc_unh_bound, aes(long, lat, group=group), 
                 fill=NA, colour="red") +
    theme(plot.background=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          plot.margin=unit(c(0,0,0,0), "in"))
  )
property_map <- ggplot() + ms_fonts +
  geom_polygon(data=lc_unh_bound, aes(long, lat, group=group), 
               fill="gray97", colour=NA) +
  geom_sf(data=towns, fill=NA, colour="gray80", size=0.2) +
  geom_polygon(data=lc_unh_bound, aes(long, lat, group=group), 
               fill=NA, colour="red", size=1) +
  geom_sf(data=lc_unh.sf, colour=NA) + # just coerces coordinates
  geom_tile(data=filter(lc_unh.df, !is.na(Property)),
            aes(lon, lat, fill=Property)) +
  scale_fill_brewer("Managed property", type="qual", palette="Paired") +
  scale_x_continuous(breaks=c(-71.4, -71.2, -71, -70.8)) +
  labs(x="Longitude", y="Latitude") +
  theme(legend.key.size=unit(0.12, "in"),
        legend.spacing=unit(0, "in")) +
  annotation_custom(grob=context_map, 
                    xmin=791464, xmax=805500,
                    ymin=4763575, ymax=4782500) +
  annotation_scale(location="bl", pad_x=unit(1.01, "in"), pad_y=unit(0.18, "in"),
                   height=unit(0.05, "in"), width_hint=0.09, 
                   bar_cols=c("gray30", "gray90"), text_cex=0.6) 
ggsave(paste0(ms_fig.dir, "UNH_map.jpeg"), property_map,
       width=140, height=70, dpi=500, units="mm")




########
## mod_output
########
N.0 <- readRDS(paste0("data/inits/N_2018_", res, ".rds")) 
lab.coord <- c(x=703000, y=4875000)
scales::trans_new("log_p1", 
                  function(x) {log(x+1)},
                  function(x) {exp(x)-1})
full.2018 <- lc.df %>% 
  mutate(N.adult=rowSums(N.0[,,7]),
         N.juv=apply(N.0[,,-7], 1, sum, na.rm=T),
         B=readRDS(paste0("data/inits/B_2018_", res, ".rds")),
         nSeed=readRDS(paste0("data/inits/nSd_2018_", res, ".rds")),
         nSdStay=readRDS(paste0("data/inits/nSdStay_2018_", res, ".rds")),
         D=readRDS(paste0("data/inits/D_2018_", res, ".rds")),
         # nFl=readRDS(paste0("data/inits/nFl_2018_", res, ".rds")),
         p.est=c(as.matrix(lc.df[,4:9]) %*% gbPopMod::set_g_p()$p))
base_map <- ggplot(full.2018, aes(lon, lat)) +
  theme_bw() + ms_fonts + labs(x="", y="") + xlim(699856, 900000) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        legend.background=element_blank(),
        legend.position=c(.87,.22),
        legend.key.height=unit(0.1, "in"),
        legend.key.width=unit(0.09, "in"))
p.N_ad <- base_map + geom_tile(aes(fill=N.adult)) + 
  scale_fill_viridis("", option="B", limits=c(0,NA), trans="log1p",
                     breaks=c(0, exp(log(max(full.2018$N.adult))/2), 
                              max(full.2018$N.adult)),
                     labels=c(0, "1.6e2", "2.6e4")) +
  ggtitle("Adult abundance") + 
  annotate("text", label="a.", x=lab.coord[1], y=lab.coord[2])
p.p_est <- base_map + geom_tile(aes(fill=p.est)) + 
  scale_fill_viridis("", option="B", limits=c(0,0.5), breaks=c(0,0.25,0.5)) + 
  ggtitle("Seedling establishment rate") + 
  annotate("text", label="b.", x=lab.coord[1], y=lab.coord[2])
p.propagule <- base_map + geom_tile(aes(fill=D+nSdStay)) + 
  scale_fill_viridis("", option="B", limits=c(0,NA), trans="log1p",
                     breaks=c(0, 
                              exp(log(max(full.2018$D + full.2018$nSdStay))/2),
                              max(full.2018$D + full.2018$nSdStay)),
                     labels=c(0, "6.8e4", "4.6e7")) + 
  ggtitle("Total propagule pressure") + 
  annotate("text", label="c.", x=lab.coord[1], y=lab.coord[2])
p.pcap_seed <- base_map + geom_tile(aes(fill=nSeed/N.adult)) + 
  scale_fill_viridis("", option="B", limits=c(0,NA), trans="log1p",
                     breaks=c(0, 
                              exp(log(max(full.2018$nSeed/(full.2018$N.adult+1)))/2),
                              max(full.2018$nSeed/(full.2018$N.adult+1))),
                     labels=c(0, "4.6e1", "2.1e3")) + 
  ggtitle("Per capita seed production") + 
  annotate("text", label="d.", x=lab.coord[1], y=lab.coord[2])
p.immigr <- base_map + geom_tile(aes(fill=D)) + 
  scale_fill_viridis("", option="B", limits=c(0,NA), trans="log1p",
                     breaks=c(0, exp(log(max(full.2018$D))/2), max(full.2018$D)),
                     labels=c(0, "1.7e4", "3.0e6")) + 
  ggtitle("Number of immigrant seeds") + 
  annotate("text", label="e.", x=lab.coord[1], y=lab.coord[2])
p.pr_immigr <- base_map + geom_tile(aes(fill=D/(nSdStay+D))) + 
  scale_fill_viridis("", option="B", limits=c(0,1), breaks=c(0, 0.5, 1)) + 
  ggtitle("Proportion of immigrant seeds") + 
  annotate("text", label="f.", x=lab.coord[1], y=lab.coord[2])
mod_output <- grid.arrange(p.N_ad, p.p_est, p.propagule, 
                           p.pcap_seed, p.immigr, p.pr_immigr, 
                           nrow=2)
ggsave(paste0(ms_fig.dir, "model_output.jpeg"), mod_output, 
       width=190, height=125, dpi=700, units="mm")
       # width=7.5, height=5, dpi=300, units="in")




########
## GSA_relInf
########
ri.df <- read_csv(paste0("out/", res, "/BRT_RI.csv")) %>%
  mutate(response=as.factor(response))
resp_col <- c(pOcc="#005a32", pSB="#74c476", pK="#99000d",
              meanNg0="#084594", medNg0="#6baed6", sdNg0="#4a1486")
ri.df$resp_pretty <- factor(ri.df$response,
                            levels=names(resp_col),
                            labels=c("Cells occupied\n(adults)", 
                                     "Cells occupied\n(seed bank)", 
                                     "Occupied cells\nat K", 
                                     "Mean\nabundance", 
                                     "Median\nabundance",
                                     "Variance\nin abundance"))
ri.df$response <- lvls_reorder(ri.df$response, 
                               match(names(resp_col), levels(ri.df$response)))
ri.df$param_pretty <- factor(ri.df$param,
                             labels=c("eta", "g[B]", "gamma", "K", "m",
                                      "mu", "N[0]", "n[ldd]", "p", "c", "f", 
                                      "s[B]", "s[c]", "s[M]", "s[N]", 
                                      "sdd[max]", "r"))
relInf <- ggplot(filter(ri.df, smp==max(ri.df$smp) & td==max(ri.df$td)), 
       aes(x=param_pretty, y=rel.inf, fill=response)) + 
  ms_fonts + coord_flip() + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") + 
  scale_fill_manual(values=resp_col) + 
  geom_bar(stat="identity", colour="gray30") + 
  facet_wrap(~resp_pretty, ncol=6) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1")) + 
  scale_x_discrete(breaks=levels(ri.df$param_pretty), 
                   labels=parse(text=levels(ri.df$param_pretty))) +
  labs(y="Relative Influence", x="Parameter") +
  theme(legend.position="none", 
        strip.text=element_text(size=9),
        panel.grid.minor=element_blank())
ggsave(paste0(ms_fig.dir, "GSA_relInf.jpeg"), relInf, 
       width=190, height=87, dpi=500, units="mm")
       # width=6.5, height=3, dpi=300, units="in")




########
## mgmt_N_B
########
out.property <- read_csv(paste0("mgmt/out/", res, "_out_property.csv")) %>%
  select(mgmt, year, Property, N_adult, B, nCell) %>%
  filter(year != 1) %>%
  mutate(year=year-1) %>%
  rename(Adults=N_adult, `Seed bank`=B) %>%
  gather(Stage, Abundance, 4:5) 
out.property$mgmt <- factor(out.property$mgmt,
                            levels=c("none", "stated", "reality", "aggressive"),
                            labels=c("No action", "Stated",
                                     "Actual", "Aggressive"))
out.compare <- out.property %>% 
  mutate(N.none=c(rep(filter(out.property, 
                             mgmt=="No action" & Stage=="Adults")$Abundance, 4),
                  rep(filter(out.property, 
                             mgmt=="No action" & Stage=="Seed bank")$Abundance, 4)))
mgmt_N_B <- ggplot(out.property) + geom_hline(yintercept=0, size=0.1) + 
  geom_line(aes(year, Abundance/nCell, group=Property, colour=Property)) + 
  facet_grid(Stage~mgmt, scales="free_y") + ms_fonts +
  scale_colour_brewer("Managed property", type="qual", palette="Paired") +
  labs(y=paste("Mean density per cell")) +
  theme(panel.background=element_rect(fill="gray95"), 
        panel.grid=element_blank(),
        legend.key.size=unit(0.12, "in"))
mgmt_N_B_diff <- ggplot(filter(out.compare, mgmt != "No action"), 
                        aes(x=year, y=Abundance/N.none, 
                            colour=Property, group=Property)) + 
  geom_line() + facet_grid(Stage~mgmt) + ms_fonts +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0.0", "0.5", "1.0")) +
  scale_colour_brewer("Managed property", type="qual", palette="Paired") +
  labs(y="Proportion of un-managed abundance") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray95"),
        legend.key.size=unit(0.12, "in"))
ggsave(paste0(ms_fig.dir, "mgmt_N_B.jpeg"), mgmt_N_B,
       width=190, height=87, dpi=500, units="mm")
       # width=6.5, height=3, dpi=300, units="in")
ggsave(paste0(ms_fig.dir, "mgmt_N_B_diff.jpeg"), mgmt_N_B_diff,
       width=190, height=87, dpi=500, units="mm")
       # width=6.5, height=3, dpi=300, units="in")






########
## Specific calculations
########
out.property <- read.csv("mgmt/out/20ac_out_property.csv")
mgmt_comp.df <- read.csv("mgmt/out/20ac_mgmt_comp.csv")

# minimum density of buckthorn
out.property %>% filter(mgmt=="aggressive") %>% group_by(Property) %>%
  summarise(min=min(N_adult/nCell)) %>% summary 

