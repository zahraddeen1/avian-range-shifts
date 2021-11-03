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
  select(aou, species_code, climate_vol, ssi, shannonE_diet, mean_area, mean_occ, Trend, logMass,migclass,
         log_Brange_Area, overlap) %>%
  filter(species_code != "wesblu") %>%
  na.omit()

trait_pca <- prcomp(all_vars[, c(3:5,9:11)], center = T, scale = T)
summary(trait_pca)

ggbiplot(trait_pca, labels = all_vars$species_code) + theme_classic(base_size = 12)
ggsave("figures/trait_biplot.pdf", units = "in", height = 8, width= 10)

## Pairwise plots of variables in SEM

# diet x hab, clim x hab
# occ x area (don't need pop trend since it is other data product?)

theme_set(theme_classic(base_size = 15))
diet_hab <- ggplot(all_vars, aes(x = shannonE_diet, y = ssi)) + 
  geom_point() +
  geom_text(data = filter(all_vars, species_code %in% c("amecro", "seaspa")), 
            aes(y = ssi - 0.1, label = species_code)) +
  xlim(-0.1, 0.9) +
  annotate(geom = "text", x = 0.8, y = 0.5, label = "Generalist") +
  annotate(geom = "text", x = 0, y = 3, label = "Specialist") +
  annotate(geom = "text", x = 0.75, y = 3, 
           label = paste0("r = ", round(cor(all_vars$ssi, all_vars$shannonE_diet, use = "pairwise.complete.obs"), 2))) +
  labs(x = "Diet niche breadth", y = "Habitat specialization")

clim_hab <- ggplot(all_vars, aes(x = climate_vol, y = ssi)) + 
  geom_point() +
  xlim(-0.05, 3.3) +
  geom_text(data = filter(all_vars, species_code %in% c("larbun", "bushti")),
            aes(y = ssi - 0.1, label = species_code)) +
  annotate(geom = "text", x = 3, y = 0.5, label = "Generalist") +
  annotate(geom = "text", x = 0.2, y = 3, label = "Specialist") +
  annotate(geom = "text", x = 3, y = 3, 
           label = paste0("r = ", round(cor(all_vars$ssi, all_vars$climate_vol, use = "pairwise.complete.obs"), 2))) +
  labs(x = "Climate niche breadth", y = "")

occ_area <- ggplot(filter(all_vars, species_code != "bushti"), aes(x = mean_occ, y = mean_area)) + 
  geom_point() +
  geom_text(data = filter(all_vars, species_code %in% c("fiscro", "cacwre")),
            aes(y = mean_area - 0.1, label = species_code)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-0.3, 0.4) +
 annotate(geom = "text", x = 0.35, y = 3.5, 
           label = paste0("r = ", round(cor(all_vars$mean_occ, all_vars$mean_area, use = "pairwise.complete.obs"), 2))) +
  labs(x = expression(paste(Delta, "Range occupancy")), y = expression(paste(Delta, "Range area")))

cowplot::plot_grid(diet_hab, clim_hab, occ_area, nrow = 2, labels = c("a", "b", "c"))
ggsave("figures/model_input_correlations.pdf", units = "in", height = 7, width = 9)

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

## Variance partitioning by order, family, genus, species

all_taxo <- clim %>%
  filter(aou %in% all_vars$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  select(aou, IOCOrder, family, genus, species, phylo_name)

vars_all_phylo <- all_vars %>%
  left_join(all_taxo)

library(nlme)

hab_nest <- lme(ssi ~ 1, random = ~ 1|family/genus, data = vars_all_phylo)
varcomp(hab_nest, T, F)

clim_nest <- lme(climate_vol ~ 1, random = ~ 1|family/genus, data = vars_all_phylo)
varcomp(clim_nest, T, F)

diet_nest <- lme(shannonE_diet ~ 1, random = ~ 1|family/genus, data = vars_all_phylo)
varcomp(diet_nest, T, F)

## Make line graph of this for supplement



## Variance partitioning: niche vs migclass/body size

varpart <- all_vars %>%
  summarize(trend_all = summary(lm(Trend ~ climate_vol + ssi + shannonE_diet + migclass + logMass, .))$r.squared,
            trend_niche = summary(lm(Trend ~ climate_vol + ssi + shannonE_diet, .))$r.squared,
            trend_other = summary(lm(Trend ~ migclass + logMass, .))$r.squared,
            occ_all = summary(lm(mean_occ ~ climate_vol + ssi + shannonE_diet + migclass + logMass, .))$r.squared,
            occ_niche = summary(lm(mean_occ ~ climate_vol + ssi + shannonE_diet, .))$r.squared,
            occ_other = summary(lm(mean_occ ~ migclass + logMass, .))$r.squared,
            area_all = summary(lm(mean_area ~ climate_vol + ssi + shannonE_diet + migclass + logMass, .))$r.squared,
            area_niche = summary(lm(mean_area ~ climate_vol + ssi + shannonE_diet, .))$r.squared,
            area_other = summary(lm(mean_area ~ migclass + logMass, .))$r.squared) %>%
  mutate(trend_niche_only = trend_all - trend_other,
         trend_other_only = trend_all - trend_niche,
         trend_shared = trend_niche_only + trend_other_only - trend_all,
         area_niche_only = area_all - area_other,
         area_other_only = area_all - area_niche,
         area_shared = area_all - area_niche_only - area_other_only,
         occ_niche_only = occ_all - occ_other,
         occ_other_only = occ_all - occ_niche,
         occ_shared = occ_all - occ_niche_only - occ_other_only) 

varpart_plot <- varpart %>%
  select(trend_niche_only:occ_shared) %>%
  pivot_longer(values_to = "r2", names_to = "var", trend_niche_only:occ_shared) %>%
  mutate(response = word(var, 1, 1, sep = "_"),
         preds = word(var, 2, sep = "_"))

theme_set(theme_classic(base_size = 15))
ggplot(varpart_plot, aes(x = response, y = r2, fill = preds)) + geom_bar(position = "stack", stat = "identity") +
  coord_flip() +
  labs(x = "Variance explained", y = "Response variable", fill = "") +
  scale_fill_viridis_d(labels = c("niche" = "Niche", "other" = "Life history", "shared" = "Shared")) +
  scale_x_discrete(labels = c("trend" = "Population trend", "occ" = expression(paste(Delta, "Range occupancy")),
                              "area" = expression(paste(Delta, "Range area")))) +
  theme(legend.position = c(0.8, 0.6))
ggsave("figures/variance_paritioning.pdf")

