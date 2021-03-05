### Get breeding range land cover composition
### Run on longleaf HPC as job array, 50 spp per job

library(sf)
library(raster)
library(dplyr)
library(stringr)

# Array ID
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# top level directory
dir <- "/proj/hurlbertlab/"

# breeding range data
range_dir <- "/proj/hurlbertlab/bird_range_shps/All/"

# land cover data
lc_dir <- "/proj/hurlbertlab/cec_north_america/north_america_2015/"

# output directory
output_dir <- "/proj/hurlbertlab/gdicecco/ch3_hab_niche/breedrange_lc/"

### Breeding range species list

# Match BBS taxonomy with breeding range polygon taxonomy
fix_mismatch <- read.csv("fix_breedingrange_genus_mismatch.csv", stringsAsFactors = F) %>%
  dplyr::select(-file, -spp_name) %>%
  filter(!is.na(old_genus)) %>%
  mutate(new_binomial = paste(old_genus, species))

bbs_spp <- read.csv("species_list.csv", stringsAsFactors = F) %>%
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

if(task_id < 8) {
  start <- as.numeric(task_id)*50 - 49
  end <- as.numeric(task_id)*50
} else {
  start <- as.numeric(task_id)*50 - 49
  end <- nrow(spp_list)
}

spp_list_subs <- spp_list[start:end, ]

### Land cover data

na_lc <- raster(paste0(lc_dir, "NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif"))

lc_crs <- st_crs(na_lc)

### For each species, save proportion cover for each land cover class across range

for(s in spp_list_subs$file) {
  
  spp <- spp_list_subs$file_binomial[spp_list_subs$file == s]
  
  br <- read_sf(paste0(range_dir, s))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
    st_transform(lc_crs)
  
  lc_crop <- crop(na_lc, breeding_range)
  lc_mask <- mask(lc_crop, breeding_range)
  
  lc_freq <- freq(lc_mask)
  lc_df <- data.frame(lc_freq)
  lc_df$spp <- s
  
  write.csv(lc_df, paste0(output_dir, spp, "_breedingrange_lc.csv"), row.names = F)
  
}
