# This script details the known occurrences of glossy buckthorn, including
# historical records from herbaria and more recent records from EDDMapS.


########
## Set up
########
make_maps <- FALSE
library(tidyverse); library(gbPopMod); library(sf); library(maps)

# land cover grid: study extent
res <- c("20ac", "9km2")[1]
load(paste0("data/USDA_", res, ".rda"))
grd.border <- lc.df %>% filter(inbd) %>%
  select(lon, lat, inbd) %>%
  raster::rasterFromXYZ() %>% 
  raster::rasterToPolygons(dissolve=T) %>%
  fortify %>% select(long, lat, group) %>%
  st_as_sf(coords=c("long", "lat"))
st_crs(grd.border) <- 32618
grd.border <- grd.border %>% st_transform(crs=4326)

# new hampshire & maine boundaries
nhme <- map("county", plot=F, fill=T) %>% st_as_sf() %>% 
  filter(grepl("maine|new hampshire",  ID))
states <- map("state", plot=F, fill=T) %>% st_as_sf()
st_crs(nhme) <- st_crs(states) <- 4326
nhme.grd <- st_crop(nhme, grd.border)

# load records
eddmaps <- read.csv("data/gb/eddmaps_gb.csv") %>%
  filter(!is.na(Longitude_Decimal) & !is.na(Latitude_Decimal)) %>%
  st_as_sf(coords=c("Longitude_Decimal", "Latitude_Decimal"))
eddmaps.nhme <- read_csv("data/gb/eddmaps_ME.csv") %>%
  select(ObsDate, Latitude, Longitude) %>%
  bind_rows(read_csv("data/gb/eddmaps_NH.csv") %>%
              select(ObsDate, Latitude, Longitude)) %>%
  filter(complete.cases(.)) %>%
  st_as_sf(coords=c("Longitude", "Latitude")) %>%
  mutate(Year=lubridate::year(as.Date(ObsDate, "%m/%d/%Y")))
herb.record <- read.csv("data/gb/gb_herbariumRecords.csv") %>%
  filter(!is.na(Longitude) & !is.na(Latitude)) %>%
  st_as_sf(coords=c("Longitude", "Latitude"))
st_crs(herb.record) <- st_crs(eddmaps) <- st_crs(eddmaps.nhme) <- 4326

# crop to study extent
eddmaps.nhme <- st_crop(eddmaps.nhme, grd.border)
herb.nhme <-  st_crop(herb.record, grd.border)

# combine EDDMapS + herbarium records
spread.nhme <- herb.nhme %>% 
  select(CollectionYear) %>% 
  rename(Year=CollectionYear) %>%
  rbind(., select(eddmaps.nhme, Year)) %>%
  # st_transform(crs=32618) %>%
  arrange(Year) 
spread.coords <- st_coordinates(st_transform(spread.nhme, crs=32618))
spread.pts <- sapply(1:nrow(spread.nhme),
                    function(x) get_pt_id(lc.df, spread.coords[x,]))
out_bds <- which(map_int(spread.pts, length)==0)
spread.pts <- data.frame(id=unlist(spread.pts), 
                         Year=spread.nhme$Year[-out_bds]) %>%
  group_by(id) %>% summarise(MinObsYear=min(Year))
lc.df$MinObsYear <- NA
lc.df$MinObsYear[spread.pts$id] <- spread.pts$MinObsYear
write_csv(lc.df, paste0("data/gb/spread_", res, ".csv"))




########
## Basic exploration
########
# Herbarium record summary: NH + ME
range(herb.nhme$CollectionYear)
sort(herb.nhme$CollectionYear)
herb.nhme %>% arrange(CollectionYear) %>% print.AsIs
plot(sort(herb.nhme$CollectionYear), 1:nrow(herb.nhme))

# First record in study extent
N.1922 <- herb.nhme %>% filter(CollectionYear==1922) %>% 
  st_transform(crs=32618) %>% st_coordinates()

# Convex hulls of spread
Years <- unique(spread.nhme$Year)[-(1:2)]
spread.hull <- vector("list", length(Years))
for(i in seq_along(Years)) {
  spread.hull[[i]] <- spread.nhme %>%
    filter(Year <= Years[i]) %>%
    st_union() %>%
    st_convex_hull()
}

spread.sf <- do.call("rbind", purrr::map(spread.hull, st_sf)) %>%
  mutate(Year=Years) %>%
  rename(geometry=.x..i..)





########
## Identify occupancy dates by pixel
########
lc.sf <- lc.df %>% filter(inbd) %>%
  st_as_sf(coords=c("lon", "lat")) %>%
  st_set_crs(32618)
hull.intersect <- st_intersects(lc.sf, st_transform(spread.sf, crs=32618)) %>%
  purrr::map_dbl(., min)
hull.intersect[is.infinite(hull.intersect)] <- NA
lc.df <- mutate(lc.df, MinObsYear=NA)
lc.df$MinObsYear[lc.df$inbd] <- Years[hull.intersect]

write_csv(lc.df, paste0("data/gb/spread_", res, ".csv"))





########
## visualization
########
if(make_maps) {
  library(viridis)
  
  # Map: National
  ggplot() + scale_colour_viridis(option="C") +
    geom_sf(data=states, fill="gray30", colour=1) +
    geom_sf(data=eddmaps, alpha=0.5) +
    geom_sf(data=herb.record, aes(colour=CollectionYear))
  
  # Map: NH + ME 
  ggplot() + scale_colour_viridis(option="C") +
    geom_sf(data=nhme.grd, fill="gray30", colour=1) +
    geom_sf(data=eddmaps.nhme) +
    geom_sf(data=herb.nhme, aes(colour=CollectionYear), size=3)
  
  # Map: NH + ME all records
  ggplot() + scale_colour_viridis(option="C") +
    geom_sf(data=nhme.grd, fill="gray30", colour=1) +
    geom_sf(data=spread.nhme, aes(colour=Year))
  
  # Map of hulls
  ggplot() + 
    geom_sf(data=nhme.grd, fill="gray30", colour=1) +
    geom_sf(data=spread.sf, fill="white", colour="white", alpha=0.05)
  
  # Map of rasterized earliest inferred date
  ggplot(lc.df) + geom_tile(aes(lon, lat, fill=MinObsYear)) + 
    geom_sf(data=st_transform(spread.nhme, crs=32618), aes(colour=Year)) +
    scale_colour_gradient(low="darkred", high="white") +
    scale_fill_gradient(low="darkred", high="white")
}










