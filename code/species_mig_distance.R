## Calculate migratory distance using range maps

library(sf)
library(tidyverse)

## Set up

info <- Sys.info()

biodrive <- ifelse(info[["sysname"]] == "Windows", "\\\\ad.unc.edu\\bio\\HurlbertLab\\", "/Volumes/HurlbertLab/")

## species range maps directory
range_dir <- paste0(biodrive, "GIS/birds/All/All/")

## Species

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

## Mig distance
## Centroid of breeding/resident range
## Centroid of wintering/resident range
## Great circle distance

mig_dist <- spp_list %>%
  mutate(mig_dist = purrr::map(file, ~{
    f <- .
    
    print(f)
    
    br <- read_sf(paste0(range_dir, f))
    names(br)[1:14] <- toupper(names(br)[1:14])
    
    # Use extant (PRESENCE 1-3), resident and breeding ranges (SEASONAL 1-2) only
    # Re-project to land cover CRS
    breeding_range <- br %>%
      filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
      st_combine(.)
    
    # Wintering range: SEASONAL (1, 3)
    wintering_range <- br %>%
      filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1, 3)) %>%
      st_combine(.)
    
    # Centroids of wintering and breeding ranges
    br_cent <- st_centroid(breeding_range)
    wint_cent <- st_centroid(wintering_range)
    
    st_distance(br_cent, wint_cent)
  }))

mig_unnest <- mig_dist %>%
  unnest(cols = c("mig_dist")) %>%
  mutate(mig_dist_m = as.numeric(mig_dist)) %>%
  dplyr::select(-mig_dist)
write.csv(mig_unnest, "derived_data/migratory_distance.csv", row.names = F)
