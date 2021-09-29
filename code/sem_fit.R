## Niche breadth SEM

library(tidyverse)
library(piecewiseSEM)

## Read in data

clim <- read_csv("derived_data/climate_niche_breadth.csv")
hab <- read_csv("derived_data/habitat_niche_ssi.csv")
diet <- read_csv("derived_data/diet_niche_breadth.csv")

range <- read_csv("derived_data/range_metrics_sampled.csv") %>%
  group_by(aou) %>%
  summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T))

poptrend <- read_csv("raw_data/BBS_1966-2017_core_trend_revised_v2.csv", 
                     col_types = cols(AOU = col_double())) %>%
  filter(Region == "SU1")

## Model table

range_mod <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range) %>%
  left_join(select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  filter(!is.na(mean_area) & !is.na(ssi))

## Fit piecewiseSEM

spp_data <- as.data.frame(range_mod)
# write.csv(spp_data, "derived_data/sem_mod_input.csv", row.names = F)

spp_psem <- psem(
  lm(mean_occ ~ shannonE_diet + climate_vol + ssi + Trend, data = spp_data),
  lm(Trend ~ ssi + climate_vol + shannonE_diet, data = spp_data),
  lm(mean_area ~ ssi + climate_vol + shannonE_diet + Trend, data = spp_data),
  data = spp_data)

summary(spp_psem)

spp_psem_simple <- psem(
  lm(mean_occ ~ shannonE_diet + climate_vol + ssi, data = spp_data),
  lm(Trend ~ ssi + climate_vol + shannonE_diet, data = spp_data),
  lm(mean_area ~ ssi + climate_vol + shannonE_diet, data = spp_data),
  data = spp_data)

summary(spp_psem_simple)


