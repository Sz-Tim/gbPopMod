# This script rasterizes the vector boundaries of each UNH woodlands property
# and identifies the landscape cells within each property

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
bbox_2sdd <- st_bbox(st_buffer(all_boundaries, dist=set_g_p()$sdd.max*2*285))
x_rng <- range(filter(lc.df, lon<=bbox_2sdd$xmax & lon>=bbox_2sdd$xmin)$x)
y_rng <- range(filter(lc.df, lat<=bbox_2sdd$ymax & lat>=bbox_2sdd$ymin)$y)
lc.df_2sdd <- filter(lc.df, x>=x_rng[1] & x<=x_rng[2] & y>=y_rng[1] & y<=y_rng[2])
lc.st <- raster::rasterFromXYZ(lc.df_2sdd[lc.df_2sdd$inbd,c(10,11,4:9,12:14)])
raster::crs(lc.st) <- sp::CRS('+init=EPSG:32618')
lc.st <- rasterToPolygons(lc.st) %>% st_as_sf()
lc.UNH.overlap <- st_intersects(all_boundaries, lc.st) %>%
  map(., ~lc.st$id.in[.]) %>% 
  setNames(unh.i$Property[match(boundary.f, unh.i$shpFileName)])
overlap.df <- data.frame(id.in=unlist(lc.UNH.overlap),
                         property=unlist(map2(lc.UNH.overlap, 
                                              names(lc.UNH.overlap), 
                                              ~rep(.y, length(.x)))))
lc.df_2sdd <- lc.df_2sdd %>% 
  mutate(in.UNH=id.in %in% overlap.df$id.in,
         Property=overlap.df$property[match(id.in, overlap.df$id.in)],
         id.full=id,
         id.in.full=id.in,
         id=row_number(),
         id.in=min_rank(na_if(inbd*id, 0)))

write_csv(lc.df_2sdd, paste0("data/USDA_", res, "_mgmt.csv"))
save(lc.df_2sdd, file=paste0("data/USDA_", res, "_mgmt.rda"))


ggplot() + 
  geom_tile(data=lc.df_2sdd, aes(lon, lat)) +
  geom_tile(data=filter(lc.df_2sdd, in.UNH), aes(lon, lat, fill=Property)) +
  scale_fill_viridis(option="D", discrete=TRUE, na.value=NA)

