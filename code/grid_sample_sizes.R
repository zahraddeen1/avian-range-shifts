### BBS sample sizes by grid resolution
### ID best early time window to max sample size

library(tidyverse)
library(sf)
library(tmap)
library(raster)

biodrive <- "\\\\ad.unc.edu\\bio\\HurlbertLab\\"

## BBS route data

routes <- read_csv(paste0(biodrive, "Databases/BBS/2017/bbs_routes_20170712.csv"))
weather <- read_csv(paste0(biodrive, "Databases/BBS/2017/bbs_weather_20170712.csv"))

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)

## NA map

na_map <- read_sf("raw_data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN")

## Consider 5 year time windows, starting year 1970-1980
## 1 degree grid
## Min threshold of 1-4 routes per grid cells, maps for these combinations

time_windows <- data.frame(start = c(1970:1980)) %>%
  mutate(end = start + 4)

## Fun: for starting year and ending year, filter BBS routes, snap to grid, plot sampled grid cells for min 1-4 routes per cell

grid_plots <- function(start_year, end_year){
  routes_subs <- RT1.routes %>%
    filter(year >= start_year, year <= end_year) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude)) %>%
    group_by(lat_cell, lon_cell) %>%
    summarize(n_routes = n_distinct(stateroute)) %>%
    st_as_sf(coords = c("lon_cell", "lat_cell")) %>%
    st_set_crs(4326)
  
  pts <- as(routes_subs, "Spatial")
  
  # Generate empty raster layer and rasterize points
 routes_raster <- raster(crs = crs(pts), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts, .)
  
  min1 <- tm_shape(na_map) + tm_polygons() + 
    tm_shape(routes_raster) + tm_raster("n_routes", palette = "YlGnBu", alpha = 0.8, breaks = c(seq(1,20, by = 2)), title = "Routes") +
    tm_layout(legend.position = c("left", "bottom"), title = "Min routes/cell = 1")
  
  # Generate rasters for 2-4 thresholds programmatically

}

