### BBS sample sizes by grid resolution
### ID best early time window to max sample size

library(tidyverse)
library(sf)
library(tmap)
library(raster)

## Bio drive paths

info <- Sys.info()

biodrive <- ifelse(info[["sysname"]] == "Windows", "\\\\ad.unc.edu\\bio\\HurlbertLab\\", "/Volumes/bio/HurlbertLab/")

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

## For routes spatial points, rasterize and plot on N. Am map for a given threshold of routes per cel
## pts = min routes per cell, sp_df = routes spatial points
routes_map <- function(min_rte, sp_df){
  df <- sp_df %>%
    filter(n_routes >= min_rte)
  
  pts <- as(df, "Spatial")
  
  # Generate empty raster layer and rasterize points
  routes_raster <- raster(crs = crs(pts), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts, .)
  
  plot1 <- tm_shape(na_map) + tm_polygons() + 
    tm_shape(routes_raster) + tm_raster("n_bins", palette = "YlGnBu", alpha = 0.8, breaks = c(seq(1,12, by = 2)), title = "Routes") +
    tm_layout(legend.position = c("left", "bottom"), title = paste0("Min routes/cell = ", min_rte))
  
  return(plot1)
}


## Fun: for starting year and ending year, filter BBS routes, snap to grid, plot sampled grid cells for min 1-4 routes per cell

grid_plots <- function(start_year, end_year){
  
  routes_subs <- RT1.routes %>%
    mutate(yr_bin = 5*floor(year/5)) %>%
    group_by(stateroute) %>%
    mutate(n_bins = n_distinct(yr_bin)) %>%
    filter(year >= start_year, year <= end_year) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude)) %>%
    group_by(lat_cell, lon_cell) %>%
    summarize(n_routes = n_distinct(stateroute), n_bins = max(n_bins)) %>%
    st_as_sf(coords = c("lon_cell", "lat_cell")) %>%
    st_set_crs(4326)
  
  # Generate rasters for thresholds from 1-4
  min1 <- routes_map(1, routes_subs) + tm_layout(main.title = paste0("Start ", start_year, ", End ", end_year))
  min2 <- routes_map(2, routes_subs)
  min3 <- routes_map(3, routes_subs)
  min4 <- routes_map(4, routes_subs)
  
  four_panel <- tmap_arrange(min1, min2, min3, min4, nrow = 2)

  return(four_panel)
}

## Make plots for 1970s:

pdf(paste0(getwd(), "/figures/grid_cell_maps.pdf"), height = 8, width = 10)
for(i in 1:11) {
  s <- time_windows[i, 1]
  e <- time_windows[i, 2]
  
  p <- grid_plots(s, e)
  
  print(p)
}
dev.off()

