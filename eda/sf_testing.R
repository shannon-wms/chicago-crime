library(sf)
library(raster)
Crime_2019 <- read_csv("data/crime-2019.csv")
newCrime <- na.omit(Crime_2019)

bounds <- st_read("data/nbd_bounds.shp")

pnts_sf <- st_as_sf(newCrime, coords = c('Longitude', 'Latitude'), crs = st_crs(bounds))


pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, bounds))
  , area = if_else(is.na(intersection), '', bounds$PRI_NEIGH[intersection])
) 

pnts
plot(pnts_sf)
