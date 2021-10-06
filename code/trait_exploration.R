## Trait data correlations

library(tidyverse)
library(ggbiplot)

## Read in data

clim <- read_csv("derived_data/climate_niche_breadth.csv")
hab <- read_csv("derived_data/habitat_niche_ssi.csv")
diet <- read_csv("derived_data/diet_niche_breadth.csv")

range <- read_csv("derived_data/range_metrics_sampled.csv") %>%
  group_by(aou) %>%
  dplyr::summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T))

poptrend <- read_csv("raw_data/BBS_1966-2017_core_trend_revised_v2.csv", 
                     col_types = cols(AOU = col_double())) %>%
  filter(Region == "SU1")

# MO Master Correlates - range size, elev heterogeneity?
ro_correlates <- read_csv("raw_data/Master_RO_Correlates_20110610.csv")

# Range overlap with BBS
range_overlap <- read_csv("derived_data/spp_bbs_range_overlap.csv") %>%
  select(aou, species_code, overlap)

## Model table

all_vars <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range) %>%
  left_join(select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  filter(!is.na(mean_area) & !is.na(ssi)) %>%
  left_join(ro_correlates, by = c("aou" = "AOU")) %>%
  left_join(range_overlap) %>%
  select(aou, species_code, climate_vol, ssi, shannonE_diet, mean_area, mean_occ, Trend, logMass,
         log_Brange_Area, overlap) %>%
  filter(species_code != "wesblu")

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = "pairwise.complete.obs"), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

pairs(all_vars[, 3:11],
      lower.panel = panel.cor,
      upper.panel = upper.panel)

## Do NMDS ordination of traits, plot spp in ordination space w/ color by range outcome

all_vars <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range) %>%
  left_join(select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  filter(!is.na(mean_area) & !is.na(ssi)) %>%
  left_join(ro_correlates, by = c("aou" = "AOU")) %>%
  left_join(range_overlap) %>%
  select(aou, species_code, climate_vol, ssi, shannonE_diet, mean_area, mean_occ, Trend, logMass,
         log_Brange_Area, overlap) %>%
  filter(species_code != "wesblu") %>%
  na.omit()

trait_pca <- prcomp(all_vars[, c(3:5,9:11)], center = T, scale = T)
summary(trait_pca)

ggbiplot(trait_pca, labels = all_vars$species_code) + theme_classic(base_size = 12)
ggsave("figures/trait_biplot.pdf", units = "in", height = 8, width= 10)

### Phylogenetic signal of niche measurements

tree_taxo <- read_csv("raw_data/BLIOCPhyloMasterTax.csv")

spp_taxo <- clim %>%
  filter(aou %in% all_vars$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  select(aou, phylo_name)

## read in 100 trees

library(ape)
library(phytools)

vars_phylo <- all_vars %>%
  left_join(spp_taxo)

bird_trees <- read.nexus("raw_data/birdtrees/output.nex")

hab_niche <- setNames(vars_phylo$ssi, vars_phylo$phylo_name)
clim_niche <- setNames(vars_phylo$climate_vol, vars_phylo$phylo_name)
diet_niche <- setNames(vars_phylo$shannonE_diet, vars_phylo$phylo_name)

hab_ls <- purrr::map_dfc(bird_trees, ~phylosig(., hab_niche, method = "lambda")$lambda)
clim_ls <- purrr::map_dfc(bird_trees, ~phylosig(., clim_niche, method = "lambda")$lambda)
diet_ls <- purrr::map_dfc(bird_trees, ~phylosig(., diet_niche, method = "lambda")$lambda)


hab_ks <- purrr::map_dfc(bird_trees, ~phylosig(., hab_niche, test = T)$K)
clim_ks <- purrr::map_dfc(bird_trees, ~phylosig(., clim_niche, test = T)$K)
diet_ks <- purrr::map_dfc(bird_trees, ~phylosig(., diet_niche, test = T)$K)

hab_df <- hab_ks %>%
  pivot_longer(names_to = "tree", values_to = "hab", 1:100)

clim_df <- clim_ks %>%
  pivot_longer(names_to = "tree", values_to = "clim", 1:100)

diet_df <- diet_ks %>%
  pivot_longer(names_to = "tree", values_to = "diet", 1:100)

hab_ldf <- hab_ls %>%
  pivot_longer(names_to = "tree", values_to = "hab", 1:100)

clim_ldf <- clim_ls %>%
  pivot_longer(names_to = "tree", values_to = "clim", 1:100)

diet_ldf <- diet_ls %>%
  pivot_longer(names_to = "tree", values_to = "diet", 1:100)

niche_phylo <- hab_df %>%
  left_join(clim_df) %>%
  left_join(diet_df) %>%
  left_join(hab_ldf, by = c("tree"),  suffix = c("_k", "_lambda")) %>%
  left_join(clim_ldf, by = c("tree"),  suffix = c("_k", "_lambda")) %>%  
  left_join(diet_ldf, by = c("tree"),  suffix = c("_k", "_lambda"))
# write.csv(niche_phylo, "derived_data/niche_phylotest_pvals.csv", row.names = F)
niche_phylo <- read_csv("derived_data/niche_phylotest_pvals.csv")

boxplot(niche_phylo[,2:7])
