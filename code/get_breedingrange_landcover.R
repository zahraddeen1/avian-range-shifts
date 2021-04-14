### Get breeding range land cover composition
### Run on longleaf HPC as job array, ~12 spp per job

library(sf)
library(raster)
library(dplyr)
library(stringr)

# Number of array jobs
jobs <- 32

# Array ID
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# top level directory
dir <- "/proj/hurlbertlab/"

# breeding range data
range_dir <- "/proj/hurlbertlab/bird_range_shps/All/"

# land cover data
lc_dir <- "/proj/hurlbertlab/cec_north_america/north_america_2015/"

# breeding range raster dir
raster_dir <- "/proj/hurlbertlab/gdicecco/breedrange_landcover_raster/"

# output directory
output_dir <- "/proj/hurlbertlab/gdicecco/ch3_hab_niche/breedrange_landcover/"

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
                                      grepl("Picoides dorsalis", .) ~ "Picoides tridactylus",
                                      grepl("Pica hudsonia", .) ~ "Pica pica",
                                      grepl("Melozone fusca", .) ~ "Melozone fuscus",
                                      grepl("Geothlypis formosa", .) ~ "Oporornis formosus",
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
  right_join(bbs_spp, by = c("file_binomial" = "matched_name")) %>%
  filter(english_common_name != "Northern Red Bishop", !is.na(file))

## Filter out species already processed

done_files <- list.files(output_dir)

spp_done <- data.frame(file = done_files) %>%
  mutate(species = word(file, 1,2, sep = "_")) %>%
  filter(!grepl("breedingrange", species))

spp_list <- spp_list %>%
  filter(!(matched_filename %in% spp_done$species))

## Species per job

n_spp <- round(nrow(spp_list)/jobs)

if(task_id < jobs) {
  start <- as.numeric(task_id)*n_spp - (n_spp - 1)
  end <- as.numeric(task_id)*n_spp
} else {
  start <- as.numeric(task_id)*n_spp - (n_spp - 1)
  end <- nrow(spp_list)
}

spp_list_subs <- spp_list[start:end, ]

### Land cover data

na_lc <- raster(paste0(lc_dir, "NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif"))

na_lc <- aggregate(na_lc, fact = 3, fun = modal)

lc_crs <- st_crs(na_lc)

### For each species, save proportion cover for each land cover class across range

for(s in spp_list_subs$file) {
  
  spp <- spp_list_subs$spp_name[spp_list_subs$file == s][1]
  
  print(paste(Sys.time(), "Starting", spp[1]))
  
  br <- read_sf(paste0(range_dir, s))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
    st_transform(lc_crs)
  
  lc_crop <- crop(na_lc, breeding_range)
  lc_mask <- mask(lc_crop, breeding_range)
  
  writeRaster(lc_mask, filename = paste0(raster_dir, spp[1], '.tif'), overwrite = T)
  
  gdalLog <- capture.output(gdalUtilities::gdalinfo(datasetname = paste0(raster_dir, spp[1], '.tif'), hist = TRUE))
  (bucxml <- as.numeric(sub('buckets.+', '', grep('buckets ', gdalLog, value = TRUE))))
  (minxml <- as.numeric(gsub('.+from | to.+', '', grep('buckets ', gdalLog, value = TRUE)) ))
  (maxxml <- as.numeric(gsub('.+to |:', '', grep('buckets ', gdalLog, value = TRUE))))
  (histxml <- as.numeric(strsplit(split = '[[:space:]]', gsub("^ |^  ", "", gdalLog[grep('buckets', gdalLog)+1]))[[1]]))
  
  labs <- seq(from = minxml, to = maxxml, length.out = bucxml)
  df <- data.frame(labs, nwlab = c(ceiling(labs[1]),
                                   round(labs[2:(bucxml-1)]),
                                   floor(labs[bucxml])), 
                   val = histxml)
  hist <- aggregate(df$val, by = list(df$nwlab), sum)
  
  print(paste(Sys.time(), spp[1], "gdalUtilities method written"))
  
  write.csv(hist, paste0(output_dir, spp[1], "_breedingrange_lc.csv"), row.names = F)
  
}