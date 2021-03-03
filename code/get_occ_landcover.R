### Match ebird occurrence records with land cover data
### Run on longleaf HPC

library(raster)
library(dplyr)
library(rgdal)
library(stringr)

# Array ID
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# ebird occ directory
occ_dir <- "/proj/hurlbertlab/gdicecco/ch3_hab_niche/occ/"

# directory to save output
output_dir <- "/proj/hurlbertlab/gdicecco/ch3_hab_niche/occ_lc/"

# land cover data
lc_dir <- "/proj/hurlbertlab/cec_north_america/north_america_2015/"

na_lc <- raster(paste0(lc_dir, "NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif"))

lc_crs <- crs(na_lc)

# Get species files based on array ID [1-6, batches of 60 spp; last batch 61]

files <- list.files(occ_dir)

if(task_id < 6) {
  start <- as.numeric(task_id)*60 - 59
  end <- as.numeric(task_id)*60
} else {
  start <- as.numeric(task_id)*60 - 59
  end <- 361
}

occ_files <- files[start:end]

# read in occ data, convert to spatial points

occ_crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

for(f in occ_files) {
  spp <- word(f, 1, sep = "_")
  
  occ <- read.csv(paste0(occ_dir, f))
  
  occ_sp <- SpatialPointsDataFrame(coords = occ[, 1:2],
                                   data = occ[, 3:6],
                                   proj4string = occ_crs)
  
  occ_transf <- spTransform(occ_sp, lc_crs)
  
  occ_lc <- extract(na_lc, occ_transf)
  
  occ_transf$lc <- occ_lc
  
  occ_data <- occ_transf@data
  
  occ_data$juneOcc <- rowMeans(occ_data[, 1:4], na.rm = T)
  
  res <- occ_data %>%
    group_by(lc) %>%
    summarize(meanOcc = mean(juneOcc, na.rm = T))
  
  write.csv(res, paste0(output_dir, spp, "_occ_by_landcover.csv"), row.names = F)
  
  print(paste0(Sys.time(), spp, " completed"))
}


