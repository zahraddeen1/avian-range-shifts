### Range shift estimates

library(tidyverse)
library(sf)
library(tmap)
library(raster)
library(concaveman)

biodrive <- "\\\\ad.unc.edu\\bio\\HurlbertLab\\"
biodrive <- "/Volumes/bio/HurlbertLab/"

## species range maps directory
range_dir <- paste0(biodrive, "GIS/birds/All/All/")

## Species

spp_list <- read_csv("derived_data/climate_niche_breadth.csv")

## BBS route data

routes <- read_csv(paste0(biodrive, "Databases/BBS/2017/bbs_routes_20170712.csv"))
weather <- read_csv(paste0(biodrive, "Databases/BBS/2017/bbs_weather_20170712.csv"))
counts <- read_csv(paste0(biodrive, "Databases/BBS/2017/bbs_counts_20170712.csv"))

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)

counts$stateroute <- counts$statenum*1000 + counts$route

# Filter BBS to rpid = 101, runtype = 1
counts_subs <- counts %>%
  filter(rpid == 101) %>%
  right_join(RT1.routes, by = c("countrynum", "statenum", "stateroute", "year"))

## Map

na_map <- read_sf("raw_data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN")

## For grid w/ BBS routes: list of species with at least 50% of breeding range in that area

# BBS route raster: start 1976, one route per cell

# 389 cells
routes_subs <- RT1.routes %>%
  mutate(yr_bin = 5*floor(year/5)) %>%
  group_by(stateroute) %>%
  mutate(n_bins = n_distinct(yr_bin)) %>%
  filter(year >= 1976) %>%
  mutate(lat_cell = round(latitude),
         lon_cell = round(longitude),
         cell_id = paste0(lon_cell, ",", lat_cell)) %>%
  group_by(lat_cell, lon_cell, cell_id) %>%
  summarize(n_routes = n_distinct(stateroute), n_bins = max(n_bins)) %>%
  st_as_sf(coords = c("lon_cell", "lat_cell")) %>%
  st_set_crs(4326) %>%
  filter(n_bins == max(n_bins))

pts <- as(routes_subs, "Spatial")

# Generate empty raster layer and rasterize points
routes_raster <- raster(crs = crs(pts), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90)))) %>%
  rasterize(pts, .)

## Species range map overlap with BBS route extent

spp_overlap <- spp_list %>%
  mutate(range_overlap = map_dbl(file, ~{
    f <- .
    
    br <- read_sf(paste0(range_dir, f))
    
    # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
    # Re-project to land cover CRS
    breeding_range <- br %>%
      filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))
    
    br <- as(breeding_range, "Spatial")
    
    range_raster <- raster(crs = crs(br), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
      rasterize(br, .)
    
    range_size <- sum(range_raster@data@values, na.rm = T)
    
    overlap <- mask(range_raster, routes_raster)
    
    range_overlap <- sum(overlap@data@values[, 1], na.rm = T)
    
    range_overlap/range_size
    
  }))
# write.csv(dplyr::select(spp_overlap, -climate_vol), "derived_data/spp_bbs_range_overlap.csv", row.names = F)

# 110 spp with > 20% range in BBS
spp <- spp_overlap %>%
  filter(range_overlap > 0.2)

## Plot species BBS occurrences
## Grid cell estimates
## Convex hull range size; concave range size

pdf("figures/spp_range_estimates.pdf", height = 8, width = 10)
for(s in spp$aou) {
  name <- spp$english_common_name[spp$aou == s]
  
  f <- spp$file[spp$aou == s]
  
  br <- read_sf(paste0(range_dir, f))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))
  
  # Spp BBS occurrences
  spp_counts <- counts_subs %>%
    filter(year >= 1976, year < 1981, aou == s) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell))
  
  # Spp grid cell occurrences
  spp_routes <- routes_subs %>%
    filter(cell_id %in% spp_counts$cell_id) %>%
    dplyr::select(-cell_id, -n_routes)

  if(nrow(spp_routes) >= 4) {
  pts <- as(spp_routes, "Spatial")
  
  routes_raster <- raster(crs = crs(pts), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts, .)

    # Spp convex hull
    
    convex <- concaveman(spp_routes, concavity = "Infinity")
    
    # Spp concave hull
    
    concave <- concaveman(spp_routes, concavity = 2)
    
    range_pol <- tm_shape(na_map) + tm_polygons() + tm_shape(breeding_range) + tm_polygons(col = "skyblue") + 
      tm_shape(routes_raster) + tm_raster(col = "n_bins", legend.show = F, palette = "Reds") + tm_layout(title = name)
    
    range_hulls <- tm_shape(na_map) + tm_polygons() + tm_shape(convex) + tm_polygons(col = "skyblue", alpha = 0.5) + 
      tm_shape(concave) + tm_polygons(col = "skyblue4", alpha = 0.5)
    
    range_maps <- tmap_arrange(range_pol, range_hulls, nrow =1)
    
    print(range_maps)
  }

}
dev.off()

