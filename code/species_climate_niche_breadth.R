## Climate niche breadth
## Hypervolume of species breeding range based on BIOCLIM variables

library(tidyverse)
library(purrr)
require(raster)
require(maps)
library(sf)
library(hypervolume)
library(rasterDT)

## BioArk directory
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

## species range maps directory
range_dir <- paste0(bioark, "/hurlbertlab/GIS/birds/All/All/")

## output directory
output_dir <- "derived_data/"

## Species list 
species_list <- read.csv("raw_data/species_list.csv", stringsAsFactors = F)

# Match BBS taxonomy with breeding range polygon taxonomy
fix_mismatch <- read.csv("derived_data/fix_breedingrange_genus_mismatch.csv", stringsAsFactors = F) %>%
  dplyr::select(-file, -spp_name) %>%
  filter(!is.na(old_genus)) %>%
  mutate(new_binomial = paste(old_genus, species))

bbs_spp <- species_list %>%
  mutate(binomial = paste(genus, species, sep = " ")) %>%
  filter(!grepl("unid.", english_common_name), !grepl("hybrid", english_common_name)) %>%
  mutate_at(c("binomial"), ~case_when(grepl("Colaptes auratus", .) ~ "Colaptes auratus",
                                      grepl("Junco hyemalis", .) ~ "Junco hyemalis",
                                      grepl("Setophaga coronata", .) ~ "Dendroica coronata",
                                      TRUE ~ .)) %>%
  left_join(fix_mismatch) %>%
  mutate(matched_name = ifelse(!is.na(old_genus), new_binomial, binomial),
         matched_filename = gsub(" ", "_", matched_name)) %>%
  filter(!is.na(matched_name))

range_files <- data.frame(file = list.files(range_dir)) %>%
  filter(grepl(".shp", file)) %>%
  mutate(spp_name = word(file, 1, 2, sep = "_"),
         file_binomial = gsub("_", " ", spp_name)) 

spp_list <- range_files %>%
  right_join(bbs_spp, by = c("file_binomial" = "matched_name"))

# WorldClim
climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
climatelayers_ss = climatelayers[[c(5, 10, 18)]]

# WorldClim crs
bio_crs <- st_crs(climatelayers)

# z-transform climate layers to make axes comparable
for (i in 1:nlayers(climatelayers_ss))
{
  climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
}

climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))

### Climate hypervolume function
### Input: species file path for range shapefile
### Output: climate hypervolume

climate_hypervolume <- function(species) {
  ## For each species: read in breeding range shapefile
  
  br <- read_sf(paste0(range_dir, species))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
    st_transform(bio_crs)
  
  ## Rasterize shapefile
  # fasterizeDT()
  
  null_rast <- raster(extent(br), res = res(climatelayers_ss_cropped), crs = bio_crs)
  
  br_rast <- fasterizeDT(breeding_range, null_rast)
  
  ## Range raster to df of coordinates
  br_coords <- rasterToPoints(br_rast)
  
  br_df <- data.frame(lon = br_coords[, 1], lat = br_coords[, 2])
  
  ## Extract clim vars
  br_clim <- extract(climatelayers_ss_cropped, br_df)
  
  br_nona <- na.omit(br_clim)
  
  if(nrow(br_nona) == 0) {
    return(NA) # Breeding range not in North America
  } else {
    ## Calculate hypervolume
    br_hyper <- hypervolume_gaussian(br_nona)
    vol <- get_volume(br_hyper)
    
    print(paste(Sys.time(), species, "complete"))
    
    return(vol)
  }
  

}

## Calculate climate hypervolume based on each species range map

spp_hypervol <- spp_list %>%
  filter(!is.na(file)) %>%
  filter(file != "Alauda_arvensis_22717415.shp", file != "Phylloscopus_borealis_22715316.shp") %>%
  mutate(climate_vol = map_dbl(file, ~climate_hypervolume(.))) %>%
  dplyr::select(file, spp_name, file_binomial, aou, english_common_name, family, genus, species, binomial, old_genus, 
                new_binomial, matched_filename, species_code, climate_vol)
  
write.csv(spp_hypervol,
            paste0(output_dir, "climate_niche_breadth.csv"),
            row.names = F)





