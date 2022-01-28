## Niche breadth SEM

library(tidyverse)
library(piecewiseSEM)
library(nlme)
library(ape)

## Read in data

clim <- read_csv("derived_data/climate_niche_breadth.csv")
hab <- read_csv("derived_data/habitat_niche_ssi.csv")
diet <- read_csv("derived_data/diet_niche_breadth.csv")

range <- read_csv("derived_data/range_metrics_sampled.csv") %>%
  group_by(aou) %>%
  summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T))

range_lo_spp <- read_csv("derived_data/range_metrics_sampled_lo_concav_fix_spp.csv") %>%
  group_by(aou) %>%
  summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T))

range_lo <- read_csv("derived_data/range_metrics_sampled_lo_concav.csv") %>%
  group_by(aou) %>%
  summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T)) %>%
  bind_rows(range_lo_spp)

range_hi <- read_csv("derived_data/range_metrics_sampled_hi_concav.csv") %>%
  group_by(aou) %>%
  summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T))


poptrend <- read_csv("raw_data/BBS_1966-2017_core_trend_revised_v2.csv", 
                     col_types = cols(AOU = col_double())) %>%
  filter(Region == "SU1")

## Write range areas csv
area_cc_var <- data.frame(concav = 2, delta_area = range$mean_area, aou = range$aou) %>%
  bind_rows(data.frame(concav = 1, delta_area = range_lo$mean_area, aou = range_lo$aou)) %>%    
  bind_rows(data.frame(concav = 3, delta_area = range_hi$mean_area, aou = range_hi$aou))
write.csv(area_cc_var, "derived_data/delta_area_concav_sensitivity.csv", row.names = F)

## Model table

range_mod <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range) %>%
  left_join(dplyr::select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  mutate_at(c("ssi"), ~.*-1) %>%
  filter(!is.na(mean_area) & !is.na(ssi))

range_mod_lo <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range_lo) %>%
  left_join(dplyr::select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  mutate_at(c("ssi"), ~.*-1) %>%
  filter(!is.na(mean_area) & !is.na(ssi))

range_mod_hi <-  clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range_hi) %>%
  left_join(dplyr::select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  mutate_at(c("ssi"), ~.*-1) %>%
  filter(!is.na(mean_area) & !is.na(ssi))

## Phylo trees

tree_taxo <- read_csv("raw_data/BLIOCPhyloMasterTax.csv")

spp_taxo <- clim %>%
  filter(aou %in% range_mod$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  dplyr::select(aou, phylo_name)

vars_phylo <- range_mod %>%
  left_join(spp_taxo)

bird_trees <- read.nexus("raw_data/birdtrees/output.nex")
tree1 <- bird_trees$tree_6755

## Fit piecewiseSEM

spp_data <- as.data.frame(vars_phylo)
# write.csv(spp_data, "derived_data/sem_mod_input.csv", row.names = F)

spp_data <- spp_data[which(spp_data$phylo_name %in% tree1$tip.label), ]
rownames(spp_data) <- spp_data$phylo_name

spp_psem <- psem(
  gls(mean_occ ~ shannonE_diet + climate_vol + ssi + Trend, 
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data),
  gls(Trend ~ ssi + climate_vol + shannonE_diet,
      correlation = corBrownian(0.5, tree1),
      na.action = na.omit,
      data = spp_data),
  gls(mean_area ~ ssi + climate_vol + shannonE_diet + Trend,
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data),
  data = spp_data)

summary(spp_psem)

plot(spp_psem,
     ns_dashed = T,
     node_attrs = list(
       shape = "rectangle", color = "black", x = 1:6, y = 1:3,
       width = 1))

spp_psem_simple <- spp_psem <- psem(
  gls(mean_occ ~ shannonE_diet + climate_vol + ssi + Trend, 
      na.action = na.omit,
      data = spp_data),
  gls(Trend ~ ssi + climate_vol + shannonE_diet,
      na.action = na.omit,
      data = spp_data),
  gls(mean_area ~ ssi + climate_vol + shannonE_diet + Trend,
      na.action = na.omit,
      data = spp_data),
  data = spp_data)

## Lo concavity

spp_taxo_lo <- clim %>%
  filter(aou %in% range_mod_lo$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  dplyr::select(aou, phylo_name)

vars_phylo_lo <- range_mod_lo %>%
  left_join(spp_taxo_lo)

spp_data_lo <- as.data.frame(vars_phylo_lo)

spp_data_lo <- spp_data_lo[which(spp_data_lo$phylo_name %in% tree1$tip.label), ]
rownames(spp_data_lo) <- spp_data_lo$phylo_name

spp_psem_lo <- psem(
  gls(mean_occ ~ shannonE_diet + climate_vol + ssi + Trend, 
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data_lo),
  gls(Trend ~ ssi + climate_vol + shannonE_diet,
      correlation = corBrownian(0.5, tree1),
      na.action = na.omit,
      data = spp_data_lo),
  gls(mean_area ~ ssi + climate_vol + shannonE_diet + Trend,
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data_lo),
  data = spp_data_lo)
summary(spp_psem_lo)

## Hi concavity

spp_taxo_hi <- clim %>%
  filter(aou %in% range_mod_hi$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  dplyr::select(aou, phylo_name)

vars_phylo_hi <- range_mod_hi %>%
  left_join(spp_taxo_hi)

spp_data_hi <- as.data.frame(vars_phylo_hi)

spp_data_hi <- spp_data_hi[which(spp_data_hi$phylo_name %in% tree1$tip.label), ]
rownames(spp_data_hi) <- spp_data_hi$phylo_name

spp_psem_hi <- psem(
  gls(mean_occ ~ shannonE_diet + climate_vol + ssi + Trend, 
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data_hi),
  gls(Trend ~ ssi + climate_vol + shannonE_diet,
      correlation = corBrownian(0.5, tree1),
      na.action = na.omit,
      data = spp_data_hi),
  gls(mean_area ~ ssi + climate_vol + shannonE_diet + Trend,
      na.action = na.omit,
      correlation = corBrownian(0.5, tree1),
      data = spp_data_hi),
  data = spp_data_hi)
summary(spp_psem_hi)


