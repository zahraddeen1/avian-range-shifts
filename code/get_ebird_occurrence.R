### Get eBird occurrence data - run on longleaf HPC

library(ebirdst)
library(raster)
library(dplyr)
library(rgdal)

# High level directory

dir <- "/proj/hurlbertlab/gdicecco/ch3_hab_niche"

### Function to get eBird median occurrence for species and save as flat csv
### Specific time window: June 2018

# Inputs:
# path = where to save occ csv
# spp_code = ebird species code (6 letters)

# Time window
start_time <- "2018-06-01"
end_time <- "2018-06-30"

get_ebird_occ <- function(path, spp_code) {
  sp_path <- ebirdst_download(species = spp_code, tifs_only = T, force = T)
  
  # Convert species raster to data frame
  
  # Load occurrence median raster
  occ <- load_raster("occurrence", path = sp_path)
  
  # Lower resolution from 3 km to 27 km
  sp_occ <- aggregate(occ, fact = 9, fun = mean)
  
  # Get and subset to dates within start/end period
  dates <- parse_raster_dates(sp_occ)
  
  target_dates <- which(dates > start_time & dates < end_time)
  
  # Subset occurrence raster to just time of interest
  occ_subs <- sp_occ[[target_dates]]
  
  # Reproject raster to degrees lat-lon
  occ_latlon <- projectRaster(occ_subs, crs = "+init=epsg:4326", method = "ngb")
  
  # Occ df - occurrence median by week
  occ_df <- rasterToPoints(occ_latlon)
  colnames(occ_df)[1:2] <- c("longitude", "latitude")
  
  # write occ data frame
  write.csv(occ_df, paste0(dir, path, spp_code, "_june_occ.csv"), row.names = F)
  
  # remove species file to save space
  unlink(sp_path, recursive = T)
}


### Species list

species <- read.csv(paste0(dir, "/species_list.csv"), stringsAsFactors = F)

ebird_bbs <- ebirdst::ebirdst_runs %>% 
  left_join(species, by = c("common_name" ="english_common_name")) %>%
  filter(!is.na(french_common_name))

### For each species, create an occurrence csv and save

p <- "/occ/"

for(s in ebird_bbs$species_code) {
  
  get_ebird_occ(path = p, spp_code = s)
  
  print(paste0(Sys.time(), " species complete: ", s))
  
}
