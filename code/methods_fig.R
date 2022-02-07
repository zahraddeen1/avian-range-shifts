## Methods fig showing range polygons and hulls

library(tidyverse)
library(sf)
library(tmap)
library(raster)
library(concaveman)
library(lwgeom)

theme_set(theme_classic(base_size = 15))

info <- Sys.info()

biodrive <- ifelse(info[["sysname"]] == "Windows", "\\\\ad.unc.edu\\bio\\HurlbertLab\\", "/Volumes/HurlbertLab/")

## species range maps directory
range_dir <- paste0(biodrive, "GIS/birds/All/All/")

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
  right_join(RT1.routes, by = c("countrynum", "statenum", "stateroute", "year")) %>%
  mutate_at(c("aou"), ~case_when(. >= 5670 & . <= 5690 ~ 5660, # Merge juncos
                                 . >= 4120 & . <= 4130 ~ 4120, # Merge flickers
                                 . >= 6550 & . <= 6560 ~ 6550, 
                                 TRUE ~ .)) %>% # Merge yellow-rumped warbler
  filter(aou != 4812) # Remove California scrub-jay (range maps not split)

## Map

na_map <- read_sf("raw_data/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN")

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

spp_list <- read_csv("derived_data/climate_niche_breadth.csv") %>%
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

## Example species: Louisiana waterthrush AOU = 6760

aou <- 6760

f <- spp_list$file[spp_list$aou == aou]

br <- read_sf(paste0(range_dir, f))
names(br)[1:14] <- toupper(names(br)[1:14])

# Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
# Re-project to land cover CRS
breeding_range <- br %>%
  filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))

# Spp BBS occurrences
spp_counts_t1 <- counts_subs %>%
  filter(year >= 1976, year < 1981, aou == 6760) %>%
  mutate(lat_cell = round(latitude),
         lon_cell = round(longitude),
         cell_id = paste0(lon_cell, ",", lat_cell))

# Spp grid cell occurrences
spp_routes_t1 <- routes_subs %>%
  filter(cell_id %in% spp_counts_t1$cell_id) %>%
  dplyr::select(-n_routes)

# Spp BBS occurrences
spp_counts_t2 <- counts_subs %>%
  filter(year >= 2013, year < 2018, aou == 6760) %>%
  mutate(lat_cell = round(latitude),
         lon_cell = round(longitude),
         cell_id = paste0(lon_cell, ",", lat_cell))

# Spp grid cell occurrences
spp_routes_t2 <- routes_subs %>%
  filter(cell_id %in% spp_counts_t2$cell_id) %>%
  dplyr::select(-n_routes)

  # Raster of occurrences time 1
  pts_t1 <- as(dplyr::select(spp_routes_t1, n_bins), "Spatial")
  
  routes_raster_t1 <- raster(crs = crs(pts_t1), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts_t1, .)
  
  # Raster of occurrences time 2
  pts_t2 <- as(dplyr::select(spp_routes_t2, n_bins), "Spatial")
  
  routes_raster_t2 <- raster(crs = crs(pts_t2), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
    rasterize(pts_t2, .)
  
  # Spp concave hull time 1
  
  concave <- spp_routes_t1 %>%
    concaveman(., concavity = 2) 
  
  concave_sub <- st_intersection(concave, st_make_valid(breeding_range))
  
  concave_nobr <- st_difference(concave, st_make_valid(breeding_range))
  
  concave_nobr_cast <- st_cast(concave_nobr, "POLYGON") %>%
    mutate(polygonID = row.names(.))
  
  # does concave area have at least 2 BBS routes in it?
  routes_sf <- RT1.routes %>%
    filter(stateroute %in% spp_counts_t1$stateroute) %>%
    st_as_sf(coords = c("longitude", "latitude")) %>%
    st_set_crs(4326) %>%
    st_intersection(concave_nobr_cast) %>%
    group_by(polygonID) %>%
    mutate(n_routes = n_distinct(stateroute)) %>%
    filter(n_routes >= 2)
  
  if(nrow(routes_sf) > 0) {
    add_shp <- concave_nobr_cast %>%
      filter(polygonID %in% routes_sf$polygonID)
    
    concave_all <- st_union(add_shp, st_buffer(concave_sub, 0))
    
    
  } else {
    concave_all <- concave_sub
  }

# Spp concave hull time 2

concave2 <- spp_routes_t2 %>%
  concaveman(., concavity = 2) %>%
  st_buffer(0)

concave_sub2 <- st_intersection(concave2, st_make_valid(breeding_range))

concave_nobr2 <- st_difference(concave2, st_make_valid(breeding_range))

concave_nobr_cast2 <- st_cast(concave_nobr2, "POLYGON") %>%
  mutate(polygonID = row.names(.))

# does concave area have at least 2 BBS routes in it?
routes_sf2 <- RT1.routes %>%
  filter(stateroute %in% spp_counts_t2$stateroute) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_intersection(concave_nobr_cast2) %>%
  group_by(polygonID) %>%
  mutate(n_routes = n_distinct(stateroute)) %>%
  filter(n_routes >= 2)

if(nrow(routes_sf2) > 0) {
  add_shp2 <- concave_nobr_cast2 %>%
    filter(polygonID %in% routes_sf2$polygonID)
  
  concave_all2 <- st_union(add_shp2, st_buffer(concave_sub2, 0))
  
  
} else {
  concave_all2 <- concave_sub2
}

## Range map with BBS occurrences

br_cells <- routes_subs %>%
  st_intersection(st_make_valid(breeding_range))

spp_routes_all <- br_cells %>%
  dplyr::select(cell_id, n_routes, n_bins, early_yrs, late_yrs) %>%
  dplyr::select(-n_routes) %>%
  rbind(., spp_routes_t1) %>%
  rbind(., spp_routes_t2) %>%
  distinct(cell_id, n_bins, early_yrs, late_yrs, geometry) %>%
  mutate(plot_legend = case_when(cell_id %in% spp_routes_t2$cell_id & cell_id %in% spp_routes_t1$cell_id ~ "Both",
                                 cell_id %in% spp_routes_t1$cell_id & !(cell_id %in% spp_routes_t2$cell_id) ~ "Early",
                                 cell_id %in% spp_routes_t2$cell_id & !(cell_id %in% spp_routes_t1$cell_id) ~ "Late",
                                 TRUE ~ "Neither"))

na_crop <- na_map %>%
  st_crop(c(ymin = 25, ymax = 48, xmin = -105, xmax = -70))

range_map <- tm_shape(na_crop) + tm_polygons() + tm_shape(breeding_range) + tm_polygons(col = "gray25") +
  tm_shape(filter(spp_routes_all, plot_legend != "Neither")) + tm_dots(col = "plot_legend", size = 0.5, title = "Time period", palette = "Dark2") +
  tm_shape(filter(spp_routes_all, plot_legend == "Neither")) + tm_symbols(col = "#E7298A", shape = 1, size = 0.4) +
  tm_add_legend(type = c("symbol"), labels = c("Not observed"), col = "#E7298A", shape = 1) +
  tm_layout(title = "a", inner.margins = c(0.05, 0.05, 0.05 ,0.05), scale = 1, legend.text.size = 1, legend.title.size = 1.5,
            main.title = "Louisiana Waterthrush")

## Concave hull with circles indicating areas retained outside the breeding range

concave_hull <- tm_shape(na_crop) + tm_polygons() + 
  tm_shape(breeding_range) + tm_borders(lwd = 3) +
  tm_shape(concave) + tm_polygons(col = "#D95F02", alpha = 0.5) +
  tm_shape(concave_all2) + tm_polygons(col = "#7570B3", alpha = 0.15) +
  tm_add_legend(col = c("#D95F02", "#7570B3"), labels = c("Early", "Late"), title = "Time period") +
  tm_layout(title = "b", inner.margins = c(0.05, 0.05, 0.05 ,0.05), scale = 1, legend.text.size = 1, legend.title.size = 1.5,
            main.title = "  ")

## Figure panels

panels <- tmap_arrange(range_map, concave_hull, nrow = 1)
tmap_save(panels, "figures/methods_fig_range.pdf", units = "in", height = 5, width = 10)

#### Supplemental figure: hi, med, lo concavity ####
# Use American crow: 4880

## Example species

aou <- 4880

f <- spp_list$file[spp_list$aou == aou]

br <- read_sf(paste0(range_dir, f))
names(br)[1:14] <- toupper(names(br)[1:14])

# Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
# Re-project to land cover CRS
breeding_range <- br %>%
  filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2))

# Spp BBS occurrences
spp_counts_t1 <- counts_subs %>%
  filter(year >= 1976, year < 1981, aou == 4880) %>%
  mutate(lat_cell = round(latitude),
         lon_cell = round(longitude),
         cell_id = paste0(lon_cell, ",", lat_cell))

# Spp grid cell occurrences
spp_routes_t1 <- routes_subs %>%
  filter(cell_id %in% spp_counts_t1$cell_id) %>%
  dplyr::select(-n_routes,-cell_id, -early_yrs, -late_yrs)

# Raster of occurrences time 1
pts_t1 <- as(spp_routes_t1, "Spatial")

routes_raster_t1 <- raster(crs = crs(pts_t1), vals = 0, resolution = c(1,1), ext = extent(c(-180, 180, -90, 90))) %>%
  rasterize(pts_t1, .)

# Spp concave hull concav=1,2,3

concave1 <- spp_routes_t1 %>%
  concaveman(., concavity = 1) %>%
  st_buffer(0)

concave2 <- spp_routes_t1 %>%
  concaveman(., concavity = 2) %>%
  st_buffer(0)

concave3 <- spp_routes_t1 %>%
  concaveman(., concavity = 3) %>%
  st_buffer(0)

routes_raster_id <- raster(routes_raster_t1, layer = 2)

na_crop <- na_map %>%
  st_crop(c(ymin = 19, ymax = 56, xmin = -142, xmax = -50))

map_cc1 <- tm_shape(na_crop) + tm_polygons() + 
  tm_shape(concave1) + tm_polygons(col = "gray45", alpha = 0.5) +
  tm_shape(routes_raster_id) + tm_raster("n_bins", alpha = 0.75, legend.show = F, palette = "Dark2") +
  tm_layout(main.title = c("(a) Concavity = 1"), title = "American crow", title.position = c("left","bottom"), 
            inner.margins = c(0.2,0,0.2,0),  outer.margins = c(0,0,0,0))

map_cc2 <- tm_shape(na_crop) + tm_polygons() + 
  tm_shape(concave2) + tm_polygons(col = "gray45", alpha = 0.5) +
  tm_shape(routes_raster_id) + tm_raster("n_bins", alpha = 0.75, legend.show = F, palette = "Dark2")+
  tm_layout(main.title = c("(b) Concavity = 2"), 
            inner.margins = c(0.2,0,0.2,0),  outer.margins = c(0,0,0,0))

map_cc3 <- tm_shape(na_crop) + tm_polygons() + 
  tm_shape(concave3) + tm_polygons(col = "gray45", alpha = 0.5) +
  tm_shape(routes_raster_id) + tm_raster("n_bins", alpha = 0.75, legend.show = F, palette = "Dark2")+
  tm_layout(main.title = c("(c) Concavity = 3"), 
          inner.margins = c(0.2,0,0.2,0), 
          outer.margins = c(0,0,0,0))

area_cc_var <- read_csv("derived_data/delta_area_concav_sensitivity.csv")
viol <- ggplot(filter(area_cc_var, aou != 7670), aes(x = as.factor(concav), y = delta_area)) +
  geom_violin(fill = "gray", draw_quantiles = c(0.5)) + theme_set(theme_classic(base_size = 15)) +
  labs(x = "Concavity", y = expression(paste(Delta, "Range area")), title = "(d)")

pdf("figures/suppl_concavity_examples.pdf")
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(map_cc1, vp = vplayout(1, 1))
print(map_cc2, vp = vplayout(1, 2))
print(map_cc3, vp = vplayout(2, 1))
print(viol, vp = vplayout(2, 2))
dev.off()


