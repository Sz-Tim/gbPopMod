Packages <- c("raster", "gbPopMod", "tidyverse", "magrittr", 
              "here", "doSNOW", "sf", "viridis")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))


# set parameters
res <- c("20ac", "9km2")[1]

# load landscape & woodland xlsx
load(paste0("data/USDA_", res, ".rda")) # loads landscape as lc.df
unh.i <- readxl::read_xlsx("data/UNH_woodlands/UNH_woodlands_summary.xlsx", 1)

# load GIS files
boundary.dir <- "data/UNH_woodlands/all boundaries"
boundary.f <- st_layers(boundary.dir)$name
UNH_Boundary_proj4string <- paste("+proj=tmerc", 
                                  "+lat_0=42.5", "+lon_0=-71.66666666666667", 
                                  "+k=0.999966667", "+x_0=300000", 
                                  "+y_0=0", "+datum=NAD83", "+units=us-ft", 
                                  "+no_defs", collapse=" ")

all_boundaries <- map(boundary.f, 
                      ~st_read(boundary.dir, .) %>% select("geometry")) %>%
  setNames(unh.i$Property[match(boundary.f, unh.i$shpFileName)]) %>%
  map(~st_set_crs(., value=UNH_Boundary_proj4string)) %>%
  map(~st_transform(., crs=32618)) %>%
  do.call("rbind", .) %>%
  mutate(property=boundary.f)
bbox_22km <- st_bbox(st_buffer(all_boundaries, dist=22000))
lc.df_22km <- filter(lc.df, lon<=bbox_22km$xmax & lon>=bbox_22km$xmin &
                       lat<=bbox_22km$ymax & lat>=bbox_22km$ymin)
lc.st <- raster::rasterFromXYZ(lc.df_22km[,c(10,11,4:9,12:14)])
raster::crs(lc.st) <- sp::CRS('+init=EPSG:32618')
lc.st <- rasterToPolygons(lc.st) %>% st_as_sf()
lc.UNH.overlap <- st_intersects(all_boundaries, lc.st) %>%
  map(., ~lc.st$id.in[.]) %>% 
  setNames(unh.i$Property[match(boundary.f, unh.i$shpFileName)])
overlap.df <- data.frame(id.in=unlist(lc.UNH.overlap),
                         property=unlist(map2(lc.UNH.overlap, 
                                              names(lc.UNH.overlap), 
                                              ~rep(.y, length(.x)))))
lc.df_22km$in.UNH <- lc.df_22km$id.in %in% overlap.df$id.in
lc.df_22km$Property <- overlap.df$property[match(lc.df_22km$id.in, overlap.df$id.in)]

write_csv(lc.df_22km, "data/USDA_UNH_mgmt.csv")
save(lc.df_22km, file="data/USDA_UNH_mgmt.rda")


ggplot() + 
  geom_tile(data=lc.df_22km, aes(lon, lat)) +
  geom_tile(data=filter(lc.df_22km, in.UNH), aes(lon, lat, fill=Property)) +
  scale_fill_viridis(option="D", discrete=TRUE, na.value=NA)

