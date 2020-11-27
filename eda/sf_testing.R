library(sf)
library(raster)
library(tidyverse)
Crime_2019 <- read_csv("data/crime-2019.csv")
newCrime <- na.omit(Crime_2019)

bounds <- st_read("data/nbd_bounds.shp")

plot(st_geometry(bounds))

pnts_sf <- st_as_sf(newCrime, coords = c('Longitude', 'Latitude'), crs = st_crs(bounds))

pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, bounds))
  , area = if_else(is.na(intersection), '', bounds$PRI_NEIGH[intersection])
)

pnts
data.frame(table(pnts$intersection))
bounds$total_crime <- data.frame(table(pnts$intersection))$Freq
plot(bounds["total_crime"])

