### Range shift estimates

library(tidyverse)
library(sf)
library(tmap)
library(raster)
library(concaveman)

theme_set(theme_classic(base_size = 15))

info <- Sys.info()

biodrive <- ifelse(info[["sysname"]] == "Windows", "\\\\ad.unc.edu\\bio\\HurlbertLab\\", "/Volumes/bio/HurlbertLab/")

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

pts <- as(dplyr::select(routes_subs, -cell_id), "Spatial")

# Generate empty raster layer and rasterize points
empty_raster <- raster(crs = crs(pts), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90)))
routes_raster <-  rasterize(pts, empty_raster)

## Species range map overlap with BBS route extent

spp_overlap <- spp_list %>%
  mutate(range_overlap = purrr::map(file, ~{
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
    
    range_missing <- range_size - range_overlap
    
    data.frame(overlap = range_overlap, range_excl = range_missing/range_size)
    
  })) %>%
  unnest(cols = c("range_overlap"))
# write.csv(dplyr::select(spp_overlap, -climate_vol), "derived_data/spp_bbs_range_overlap.csv", row.names = F)

# 110 spp with > 20% range in BBS
spp <- spp_overlap %>%
  filter(range_excl <= 0.8)

## Range area and range occupancy changes
## Sample 1 route per grid cell, 500x
## For each species+time window, calculate # grid cells occuppied, max and min lat & lon, range occ (occ cells/cells in polygon)

range_samp <- data.frame(sim = c(1:499)) %>%
  group_by(sim) %>%
  nest() %>%
  mutate(range_metrics = purrr::map(sim, ~{
    
    all_routes <- RT1.routes %>%
      mutate(yr_bin = 5*floor(year/5)) %>%
      group_by(stateroute) %>%
      mutate(n_bins = n_distinct(yr_bin)) %>%
      filter(year >= 1976) %>%
      mutate(lat_cell = round(latitude),
             lon_cell = round(longitude),
             cell_id = paste0(lon_cell, ",", lat_cell)) %>%
      filter(cell_id %in% routes_subs$cell_id, n_bins == 11) %>%
      distinct(stateroute, latitude, longitude, bcr, n_bins, lat_cell, lon_cell, cell_id)
    
    sample_routes <- all_routes %>%
      group_by(cell_id) %>%
      sample_n(1)
    
    # Spp BBS occurrences
    spp_counts_t1 <- counts_subs %>%
      filter(year >= 1976, year < 1981, aou %in% spp$aou) %>%
      mutate(lat_cell = round(latitude),
             lon_cell = round(longitude),
             cell_id = paste0(lon_cell, ",", lat_cell)) %>%
      filter(cell_id %in% sample_routes$cell_id) %>%
      group_by(aou) %>%
      nest() %>%
      mutate(metrics = purrr::map(data, ~{
        df <- .
        
        concave <- routes_subs %>%
          filter(cell_id %in% df$cell_id) %>%
          concaveman(., concavity = 2)
        
        data.frame(max_lat = max(df$lat_cell),
                   max_lon = max(df$lon_cell),
                   min_lat = min(df$lat_cell),
                   min_lon = min(df$lon_cell),
                   total_cells = n_distinct(df$cell_id),
                   area = st_area(concave))
      })) %>%
      dplyr::select(-data) %>%
      unnest(cols = c("metrics"))
    
    
    spp_counts_t2 <- counts_subs %>%
      filter(year >= 2013, year <= 2017, aou %in% spp$aou) %>%
      mutate(lat_cell = round(latitude),
             lon_cell = round(longitude),
             cell_id = paste0(lon_cell, ",", lat_cell)) %>%
      filter(cell_id %in% sample_routes$cell_id) %>%
      group_by(aou) %>%
      nest() %>%
      mutate(metrics = purrr::map(data, ~{
        df <- .
        
        concave <- routes_subs %>%
          filter(cell_id %in% df$cell_id) %>%
          concaveman(., concavity = 2)
        
        data.frame(max_lat = max(df$lat_cell),
                   max_lon = max(df$lon_cell),
                   min_lat = min(df$lat_cell),
                   min_lon = min(df$lon_cell),
                   total_cells = n_distinct(df$cell_id),
                   area = st_area(concave))
      })) %>%
      dplyr::select(-data) %>%
      unnest(cols = c("metrics"))
    
    spp_change <- spp_counts_t1 %>%
      left_join(spp_counts_t2, by = c("aou"), suffix = c("_t1", "_t2")) %>%
      left_join(spp_overlap, by = c("aou")) %>%
      mutate(delta_area = (area_t2 - area_t1)/1000000,
             delta_occ = (total_cells_t2 - total_cells_t1)/overlap)
    
    spp_change
  }))

range_metrics <- range_samp %>%
  dplyr::select(-data) %>%
  unnest(cols = c("range_metrics"))
write.csv(range_metrics, "derived_data/range_metrics_sampled.csv", row.names = F)

range_metrics_sum <- range_metrics %>%
  mutate(delta_lat = (max_lat_t2 - min_lat_t2) - (max_lat_t1 - min_lat_t1),
         delta_lon = (max_lon_t2 - min_lon_t2) - (max_lon_t1 - min_lon_t1)) %>%
  group_by(aou) %>%
  summarize(mean_delta_area = mean(delta_area, na.rm = T),
            sd_delta_area = sd(delta_area, na.rm = T),
            mean_delta_occ = mean(delta_occ, na.rm = T),
            sd_delta_occ = sd(delta_occ, na.rm = T),
            mean_delta_lat = mean(delta_lat, na.rm = T),
            sd_delta_lat = sd(delta_lat, na.rm = T),
            mean_delta_lon = mean(delta_lon, na.rm = T),
            sd_delta_lon = sd(delta_lon, na.rm = T))

area_occ_r <- cor(range_metrics_sum$mean_delta_area, range_metrics_sum$mean_delta_occ, use = "pairwise.complete.obs")

area_plot <- ggplot(range_metrics_sum, aes(x = mean_delta_area)) + geom_histogram(col = "white") + labs(x = "Change in area (km^2)")
occ_plot <- ggplot(range_metrics_sum, aes(x = mean_delta_occ)) + geom_histogram(col = "white") + labs(x = "Change in occupancy")
area_occ <- ggplot(range_metrics_sum, aes(x = mean_delta_area, y = mean_delta_occ)) + geom_point(alpha = 0.5) + 
  labs(x = "Change in area (km^2)", y = "Change in occupancy") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  annotate(x = 200, y = -1, geom = "text", label = paste0("r = ", round(area_occ_r, 2)))
cowplot::plot_grid(area_plot, occ_plot, area_occ, nrow = 2)
ggsave('figures/range_area_occ_hist.pdf', units = "in", height = 8, width = 10)

## Plot species BBS occurrences
## Grid cell estimates
## Convex hull range size; concave range size
## Mean, conf int of range area and range occ deltas

pdf("figures/spp_range_estimates.pdf", height = 8, width = 10)
for(s in spp$aou) {
  name <- spp$english_common_name[spp$aou == s]
  
  f <- spp$file[spp$aou == s]
  
  d_area <- range_metrics_sum$mean_delta_area[range_metrics_sum$aou == s]
  d_occ <- range_metrics_sum$mean_delta_occ[range_metrics_sum$aou == s]
  
  br <- read_sf(paste0(range_dir, f))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))
  
  # Spp BBS occurrences
  spp_counts_t1 <- counts_subs %>%
    filter(year >= 1976, year < 1981, aou == s) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell))
  
  # Spp grid cell occurrences
  spp_routes_t1 <- routes_subs %>%
    filter(cell_id %in% spp_counts_t1$cell_id) %>%
    dplyr::select(-cell_id, -n_routes)
  
  # Spp BBS occurrences
  spp_counts_t2 <- counts_subs %>%
    filter(year >= 2013, year < 2018, aou == s) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell))
  
  # Spp grid cell occurrences
  spp_routes_t2 <- routes_subs %>%
    filter(cell_id %in% spp_counts_t2$cell_id) %>%
    dplyr::select(-cell_id, -n_routes)

  if(nrow(spp_routes_t1) >= 4) {
  pts_t1 <- as(spp_routes_t1, "Spatial")
  
  routes_raster_t1 <- raster(crs = crs(pts_t1), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts_t1, .)

  pts_t2 <- as(spp_routes_t2, "Spatial")
  
  routes_raster_t2 <- raster(crs = crs(pts_t2), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts_t2, .)
  
    # Spp convex hull
    
    convex <- concaveman(spp_routes_t1, concavity = "Infinity")
    convex2 <- concaveman(spp_routes_t2, concavity = "Infinity")
    
    # Spp concave hull
    
    concave <- concaveman(spp_routes_t1, concavity = 2)
    concave2 <- concaveman(spp_routes_t2, concavity = 2)
    
    range_pol_1 <- tm_shape(na_map) + tm_polygons() + tm_shape(breeding_range) + tm_polygons(col = "skyblue") + 
      tm_shape(routes_raster_t1) + tm_raster(col = "n_bins", legend.show = F, palette = "Reds") + tm_layout(title = name)
    
    range_hulls_1 <- tm_shape(na_map) + tm_polygons() + tm_shape(convex) + tm_polygons(col = "skyblue", alpha = 0.5) + 
      tm_shape(concave) + tm_polygons(col = "skyblue4", alpha = 0.5)
    
    range_pol_2 <- tm_shape(na_map) + tm_polygons() + tm_shape(breeding_range) + tm_polygons(col = "skyblue") + 
      tm_shape(routes_raster_t2) + tm_raster(col = "n_bins", legend.show = F, palette = "Reds") + tm_layout(title = paste0("change area = ", round(d_area), "\nchange occ = ", round(d_occ, 2)))
    
    range_hulls_2 <- tm_shape(na_map) + tm_polygons() + tm_shape(convex2) + tm_polygons(col = "skyblue", alpha = 0.5) + 
      tm_shape(concave2) + tm_polygons(col = "skyblue4", alpha = 0.5)
    
    range_maps <- tmap_arrange(range_pol_1, range_hulls_1, range_pol_2, range_hulls_2, nrow =2)
    
    print(range_maps)
  }

}
dev.off()

