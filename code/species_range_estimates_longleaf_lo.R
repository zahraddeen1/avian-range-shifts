## Spp range estimates with lower concavity

## Range size estimates - iterations in parallel on longleaf

### Range shift estimates

library(tidyverse)
library(sf)
library(raster)
library(concaveman)
library(foreach)
library(doParallel)
library(lwgeom)

c <- 10
registerDoParallel(cores = c)

proj_path <- "/proj/hurlbertlab/"

## species range maps directory
range_dir <- paste0(proj_path, "bird_range_shps/All/")

## Species

spp_list <- read_csv("/proj/hurlbertlab/gdicecco/avian-range-shifts/derived_data/climate_niche_breadth.csv") %>%
  mutate_at(c("aou"), ~case_when(. >= 5670 & . <= 5690 ~ 5660, # Merge juncos
                                 . >= 4120 & . <= 4130 ~ 4120, # Merge flickers
                                 . >= 6550 & . <= 6560 ~ 6550,
                                 TRUE ~ .)) %>% # Merge yellow-rumped warbler
  mutate_at(c("english_common_name"), ~case_when(. == "(Myrtle Warbler) Yellow-rumped Warbler" | . == "(Audubon's Warbler) Yellow-rumped Warbler" ~ "Yellow-rumped Warbler",
                                                 grepl("Dark-eyed Junco", .) ~ "Dark-eyed Junco",
                                                 grepl("Flicker", .) ~ "Northern Flicker",
                                                 TRUE ~ .)) %>%
  dplyr::select(-species, -climate_vol) %>%
  distinct()

taxo <- read_csv("/proj/hurlbertlab/gdicecco/avian-range-shifts/raw_data/Bird_Taxonomy_20110227.csv")

## BBS route data

routes <- read_csv(paste0(proj_path, "bbs/bbs_routes_20170712.csv"))
weather <- read_csv(paste0(proj_path, "bbs/bbs_weather_20170712.csv"))
counts <- read_csv(paste0(proj_path, "bbs/bbs_counts_20170712.csv"))

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)

counts$stateroute <- counts$statenum*1000 + counts$route

# Filter BBS to rpid = 101, runtype = 1
counts_subs <- counts %>%
  filter(rpid == 101) %>%
  right_join(RT1.routes, by = c("countrynum", "statenum", "stateroute", "year")) %>%
  mutate_at(c("aou"), ~case_when(. >= 5670 & . <= 5690 ~ 5660, # Merge juncos
                                 . >= 4120 & . <= 4130 ~ 4120, # Merge flickers
                                 . >= 6550 & . <= 6560 ~ 6550, 
                                 TRUE ~ .)) %>% # Merge yellow-rumped warbler
  filter(aou != 4812) # Remove California scrub-jay (range maps not split)

## For grid w/ BBS routes: list of species with at least 50% of breeding range in that area

# BBS route raster: start 1976, one route per cell

# 383 cells
routes_subs <- RT1.routes %>%
  mutate(yr_bin = 5*floor(year/5)) %>%
  group_by(stateroute) %>%
  mutate(n_bins = n_distinct(yr_bin)) %>%
  group_by(stateroute, yr_bin) %>%
  mutate(n_surveys = n_distinct(year)) %>%
  filter(year >= 1976) %>%
  mutate(lat_cell = round(latitude),
         lon_cell = round(longitude),
         cell_id = paste0(lon_cell, ",", lat_cell)) %>%
  group_by(lat_cell, lon_cell, cell_id) %>%
  summarize(n_routes = n_distinct(stateroute), n_bins = max(n_bins), 
            early_yrs = max(n_surveys[yr_bin == 1975]), late_yrs = max(n_surveys[yr_bin == 2010])) %>%
  st_as_sf(coords = c("lon_cell", "lat_cell")) %>%
  st_set_crs(4326) %>%
  filter(n_bins == max(n_bins), early_yrs >= 2, late_yrs >= 2)

# BBS routes area

bbs_convex <- concaveman(routes_subs, concavity = "Infinity")

## Species range map overlap with BBS route extent

spp_overlap <- read_csv("/proj/hurlbertlab/gdicecco/avian-range-shifts/derived_data/spp_bbs_range_overlap.csv")

# 229 spp with at least 50% range in BBS convex hull
spp <- spp_overlap %>%
  mutate_at(c("overlap"), ~as.numeric(.)) %>%
  filter(overlap > 0.5)

## Fit concave hulls function, needs AOU and BBS counts data for that species
# Calculate range metrics

concave_area <- function(aou, df) {
  print(aou)
  
  f <- spp_list$file[spp_list$aou == aou]
  
  br <- read_sf(paste0(range_dir, f))
  names(br)[1:14] <- toupper(names(br)[1:14])
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))
  
  concave <- routes_subs %>%
    filter(cell_id %in% df$cell_id) %>%
    concaveman(., concavity = 1)
  
  concave_sub <- st_intersection(concave, breeding_range)
  
  concave_nobr <- st_difference(concave, breeding_range)
  
  concave_nobr_cast <- st_cast(concave_nobr, "POLYGON") %>%
    mutate(polygonID = row.names(.))
  
  bbs_br <- st_intersection(bbs_convex, breeding_range)
  
  # does concave area have at least 2 BBS routes in it?
  routes_sf <- routes %>%
    st_as_sf(coords = c("longitude", "latitude")) %>%
    st_set_crs(4326) %>%
    filter(stateroute %in% df$stateroute) %>%
    st_intersection(concave_nobr_cast) %>%
    group_by(polygonID) %>%
    mutate(n_routes = n_distinct(stateroute)) %>%
    filter(n_routes >= 2)
  
  if(nrow(routes_sf) > 0) {
    add_shp <- concave_nobr_cast %>%
      filter(polygonID %in% routes_sf$polygonID)
    
    concave_all <- st_union(add_shp, concave_sub)
    
    
  } else {
    concave_all <- concave_sub
  }
  
  # Check for outliers in concave algorithm
  # Area of concave hull that is not in BBS/Breeding Range intersection
  concave_all_noBbsBr <- st_difference(concave_all, bbs_br)
  
  # Area of concave hull + BBS/Breeding Range intersection
  concave_all_BbsBr <- st_union(concave_all_noBbsBr, bbs_br)
  
  range_outlier <- as.numeric(sum(st_area(concave_all_noBbsBr)))
  range_outlier_denom <- as.numeric(sum(st_area(concave_all_BbsBr)))
  
  return(data.frame(max_lat = max(df$lat_cell),
                    max_lon = max(df$lon_cell),
                    min_lat = min(df$lat_cell),
                    min_lon = min(df$lon_cell),
                    total_cells = n_distinct(df$cell_id),
                    area = sum(st_area(concave_all)),
                    outlier = range_outlier,
                    outlier_denom = range_outlier_denom))
}

possibly_concave_area <- possibly(concave_area, data.frame(max_lat = NA,
                                                           max_lon = NA,
                                                           min_lat = NA,
                                                           min_lon = NA,
                                                           total_cells = NA,
                                                           area = NA,
                                                           outlier = NA,
                                                           outlier_denom = NA))

## Range area and range occupancy changes
## Sample 1 route per grid cell, 500x
## For each species+time window, calculate # grid cells occuppied, max and min lat & lon, range occ (occ cells/cells in polygon)

# 1.62 hours per iteration
range_samp <- foreach(i=1:500) %dopar% {
  
  all_routes <- RT1.routes %>%
    mutate(yr_bin = 5*floor(year/5)) %>%
    group_by(stateroute) %>%
    mutate(n_bins = n_distinct(yr_bin)) %>%
    group_by(stateroute, yr_bin) %>%
    mutate(n_surveys = n_distinct(year)) %>%
    filter(year >= 1976) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell)) %>%
    group_by(cell_id) %>%
    mutate(early_yrs = max(n_surveys[yr_bin == 1975]), 
           late_yrs = max(n_surveys[yr_bin == 2010])) %>%
    filter(cell_id %in% routes_subs$cell_id, n_bins == 11, early_yrs >= 2, late_yrs >= 2) %>%
    distinct(stateroute, latitude, longitude, bcr, n_bins, lat_cell, lon_cell, cell_id)
  
  sample_routes <- all_routes %>%
    group_by(cell_id) %>%
    sample_n(1)
  
  # Spp BBS occurrences - must be present at a route >= 2/5 times to count as a presence
  spp_counts_t1 <- counts_subs %>%
    filter(year >= 1976, year < 1981, aou %in% spp$aou) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell)) %>%
    filter(stateroute %in% sample_routes$stateroute) %>%
    filter(cell_id %in% sample_routes$cell_id) %>%
    group_by(aou, cell_id, stateroute) %>%
    mutate(n_years = n_distinct(year)) %>%
    filter(n_years >= 2) %>%
    group_by(aou) %>%
    mutate(n_cells = n_distinct(cell_id)) %>%
    filter(n_cells >= 4) %>%
    nest() 
  
  
  res <- data.frame(aou = c(), max_lat = c(), max_lon = c(), min_lat = c(), min_lon = c(), 
                    total_cells = c(), area = c(), outlier = c())
  
  for(a in spp_counts_t1$aou) {
    df <- spp_counts_t1$data[spp_counts_t1$aou == a][[1]]
    
    res1 <- possibly_concave_area(a, df)
    
    res <- rbind(res, data.frame(aou = a, res1))
  }
  
  spp_counts_t2 <- counts_subs %>%
    filter(year >= 2013, year <= 2017, aou %in% spp$aou) %>%
    mutate(lat_cell = round(latitude),
           lon_cell = round(longitude),
           cell_id = paste0(lon_cell, ",", lat_cell)) %>%
    filter(stateroute %in% sample_routes$stateroute) %>%
    filter(cell_id %in% sample_routes$cell_id) %>%
    group_by(aou, cell_id, stateroute) %>%
    mutate(n_years = n_distinct(year)) %>%
    filter(n_years >= 2) %>%
    group_by(aou) %>%
    mutate(n_cells = n_distinct(cell_id)) %>%
    filter(n_cells >= 4) %>%
    nest() 
  
  res2 <- data.frame(aou = c(), max_lat = c(), max_lon = c(), min_lat = c(), min_lon = c(), 
                     total_cells = c(), area = c(), outlier = c())
  
  for(a in spp_counts_t2$aou) {
    df <- spp_counts_t2$data[spp_counts_t2$aou == a][[1]]
    
    res1 <- possibly_concave_area(a, df)
    
    res2 <- rbind(res2, data.frame(aou = a, res1))
  }
  
  res %>%
    left_join(res2, by = c("aou"), suffix = c("_t1", "_t2")) %>%
    left_join(spp_overlap, by = c("aou")) %>%
    mutate(delta_area = (as.numeric(area_t2) - as.numeric(area_t1))/as.numeric(area_t1),
           delta_occ = (total_cells_t2 - total_cells_t1)/overlap_cells)
  
}

range_metrics <- do.call(rbind.data.frame, range_samp)

write.csv(range_metrics, "/proj/hurlbertlab/gdicecco/avian-range-shifts/derived_data/range_metrics_sampled_lo_concav.csv", row.names = F)


