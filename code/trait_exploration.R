## Trait data correlations

library(tidyverse)
# library(ggbiplot)
library(nlme)
library(ape)
library(phytools)
library(cowplot)
library(vegan)

## Read in data

# migratory distance
mig_unnest <- read_csv("derived_data/migratory_distance.csv") %>%
  dplyr::select(aou, mig_dist_m)

# Taxonomy
tree_taxo <- read_csv("raw_data/BLIOCPhyloMasterTax.csv")

## read in 100 trees
bird_trees <- read.nexus("raw_data/birdtrees/output.nex")

# Model variables
clim <- read_csv("derived_data/climate_niche_breadth.csv")
hab <- read_csv("derived_data/habitat_niche_ssi.csv")
diet <- read_csv("derived_data/diet_niche_breadth.csv")

range <- read_csv("derived_data/range_metrics_sampled.csv") %>%
  group_by(aou) %>%
  dplyr::summarize(mean_area = mean((area_t2-area_t1)/area_t1, na.rm = T),
                   sd_area = sd(mean(area_t2 - area_t1)/area_t1, na.rm = T),
            mean_occ = mean((total_cells_t2 - total_cells_t1)/overlap_cells, na.rm = T),
            min_cells = min(c(total_cells_t1, total_cells_t2), na.rm = T))

range_all <- read_csv("derived_data/range_metrics_sampled.csv")

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
  filter(min_cells > 10) %>%
  left_join(select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  filter(!is.na(mean_area) & !is.na(ssi)) %>%
  left_join(ro_correlates, by = c("aou" = "AOU")) %>%
  left_join(range_overlap) %>%
  select(aou, species_code, climate_vol, ssi, shannonE_diet, mean_area, mean_occ, Trend, logMass,
         log_Brange_Area, Brange_Area_km2, overlap)

## Taxonomy included
spp_taxo <- clim %>%
  filter(aou %in% all_vars$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  select(aou, phylo_name)

all_taxo <- clim %>%
  filter(aou %in% all_vars$aou) %>%
  mutate("phylo_name" = case_when(aou == 7222 ~ "Troglodytes_troglodytes",
                                  aou == 6760 ~ "Seiurus_motacilla",
                                  aou ==  6410 ~ "Vermivora_pinus",
                                  aou == 5780 ~ "Aimophila_cassinii",
                                  TRUE ~ matched_filename)) %>%
  left_join(tree_taxo, by = c("phylo_name" = "TipLabel")) %>%
  dplyr::select(aou, IOCOrder, family, genus, species, phylo_name)

vars_phylo <- all_vars %>%
  left_join(spp_taxo)

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

all_vars_taxo <- all_vars %>%
  left_join(all_taxo) %>%
  dplyr::group_by(family) %>%
  mutate(n_spp = n_distinct(aou)) %>%
  mutate(family_plot = case_when(family %in% c("Hirundinidae", "Corvidae", "Parulidae", "Icteridae") ~ family,
                                 family == "Emberizidae" ~ "Passerellidae",
                                 TRUE ~ "Other"))

family_cols <- c(RColorBrewer::brewer.pal(5,"Set1"), "gray")

theme_set(theme_classic(base_size = 15))
diet_hab <- ggplot(all_vars_taxo, aes(x = shannonE_diet, y = -1*ssi, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha= 0.75) +
  annotate(geom = "text", y = -3, x = 0.1, label = "Specialist") +
  annotate(geom = "text", y = -3, x = 0.75, label = "Generalist") +
  annotate(geom = "text", y = 0, x = 0.1, label = "Generalist") +
  annotate(geom = "text", x = 0.75, y = -0.5, 
           label = paste0("r = ", round(cor(-1*all_vars_taxo$ssi, all_vars_taxo$shannonE_diet, use = "pairwise.complete.obs"), 2))) +
  labs(x = " ", y = "Habitat niche breadth", col = "Family", size = "Breeding range area (km^2)") +
  scale_color_manual(values = family_cols) +
  scale_size_continuous(range = c(3, 10))

clim_hab <- ggplot(all_vars_taxo, aes(x = climate_vol, y = -1*ssi, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha = 0.75) +
  annotate(geom = "text", y = -3, x = 2.85, label = "Generalist") +
  annotate(geom = "text", y = -3, x = 0.35, label = "Specialist") +
  annotate(geom = "text", y = -0.5, x = 0.35, label = "Generalist") +
  annotate(geom = "text", x = 3, y = -0.5, 
           label = paste0("r = ", round(cor(-1*all_vars$ssi, all_vars$climate_vol, use = "pairwise.complete.obs"), 2))) +
  labs(x = "Climate niche breadth", y = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = family_cols) +
  scale_size_continuous(range = c(3, 10))

diet_clim <- ggplot(all_vars_taxo, aes(x = shannonE_diet, y = climate_vol, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha= 0.75) +
  annotate(geom = "text", y = 0, x = 0.1, label = "Specialist") +
  annotate(geom = "text", y = 0, x = 0.75, label = "Generalist") +
  annotate(geom = "text", y = 3.2, x = 0.1, label = "Generalist") +
  annotate(geom = "text", x = 0.75, y = 3, 
           label = paste0("r = ", round(cor(all_vars_taxo$climate_vol, all_vars_taxo$shannonE_diet, use = "pairwise.complete.obs"), 2))) +
  labs(x = "Diet niche breadth", y = "Climate niche breadth", col = "Family", size = "Breeding range area") +
  scale_color_manual(values = family_cols) +
  theme(legend.position = "none")+
  scale_size_continuous(range = c(3, 10))

scatter_legend <- get_legend(diet_hab)

plot_grid(diet_hab + theme(legend.position = "none"), clim_hab, 
                    diet_clim, scatter_legend,
                   nrow = 2, labels = c("a", "b", "c"))
ggsave("figures/model_input_correlations.pdf", units = "in", height = 8, width = 10)

# 3 panels of response vars

trend_area <- ggplot(all_vars_taxo, aes(x = Trend, y = mean_area, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  annotate(geom = "text", x = 2.2, y = 3.5, 
           label = paste0("r = ", round(cor(all_vars$Trend, all_vars$mean_area, use = "pairwise.complete.obs"), 2))) +
  labs(x = "Population trend", y = expression(paste(Delta, "Range area"))) +
  theme(legend.position = "none") +
  scale_color_manual(values = family_cols)+
  scale_size_continuous(range = c(3, 8))

occ_area <- ggplot(filter(all_vars_taxo, species_code != "bushti"), aes(x = mean_occ, y = mean_area, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-0.3, 0.4) +
  annotate(geom = "text", x = 0.35, y = 3.5, 
           label = paste0("r = ", round(cor(all_vars$mean_occ, all_vars$mean_area, use = "pairwise.complete.obs"), 2))) +
  labs(x = expression(paste(Delta, "Range occupancy")), y = "") +
  scale_color_manual(values = family_cols) +
  theme(legend.position = "none")+
  scale_size_continuous(range = c(3, 8))

trend_occ <- ggplot(all_vars_taxo, aes(x = Trend, y = mean_occ, col = fct_relevel(family_plot, "Other", after = Inf), size = Brange_Area_km2)) + 
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  annotate(geom = "text", x = 2.2, y = 0.35, 
           label = paste0("r = ", round(cor(all_vars$mean_occ, all_vars$Trend, use = "pairwise.complete.obs"), 2))) +
  labs(y = expression(paste(Delta, "Range occupancy")), x = "Population trend",  col = "Family", size = "Breeding range area") +
  scale_color_manual(values = family_cols)+
  scale_size_continuous(range = c(3, 8))

response_legend <- get_legend(trend_occ)
  
plot_grid(trend_area, occ_area, trend_occ + theme(legend.position = "none"), response_legend, nrow =2,
          labels = c("a", "b", "c"))
ggsave("figures/model_response_correlations.pdf", units = "in", height= 8, width = 10)

cor.test(all_vars$Trend, all_vars$mean_area, use = "pairwise.complete.obs")
cor.test(all_vars$mean_occ, all_vars$mean_area, use = "pairwise.complete.obs")
cor.test(all_vars$mean_occ, all_vars$Trend, use = "pairwise.complete.obs")

### Phylogenetic signal of niche measurements

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

vars_all_phylo <- all_vars %>%
  left_join(all_taxo) 

library(nlme)

hab_nest <- lme(ssi ~ 1, random = ~ 1|family/genus, data = vars_all_phylo)
hab_var <- varcomp(hab_nest, T, F)

clim_nest <- lme(climate_vol ~ 1, random = ~ 1|family/genus, data = vars_all_phylo)
clim_var <- varcomp(clim_nest, T, F)

diet_nest <- lme(shannonE_diet ~ 1, random = ~ 1|family/genus, data = filter(vars_all_phylo, !is.na(shannonE_diet)))
diet_var <- varcomp(diet_nest, T, F)

## Make line graph of this for supplement

niche_var <- data.frame(niche = "hab", family = hab_var[[1]], 
                        genus =hab_var[[2]],
                        species = hab_var[[3]]) %>%
  bind_rows(data.frame(niche = "clim", family = clim_var[[1]], 
                       genus =clim_var[[2]],
                       species = clim_var[[3]])) %>%
  bind_rows(data.frame(niche = "diet", family = diet_var[[1]], 
             genus =diet_var[[2]],
             species = diet_var[[3]])) %>%
  pivot_longer(names_to = "phylo", values_to = "var", family:species)

theme_set(theme_classic(base_size = 15))
ggplot(niche_var, aes(x = phylo, y = var, col = niche, group = niche)) + 
  geom_line(cex = 1) +
  scale_x_discrete(labels = c("family" = "Family",
                              "genus" = "Genus",
                              "species" = "Species")) +
  scale_color_brewer(palette = "Dark2",
                     labels = c("clim" = "Climate",
                                "diet" = "Diet",
                                "hab" = "Habitat")) +
  labs(x = "", y = "Variance explained in niche breadth", col = "Niche axis") +
  theme(legend.position = c(0.8, 0.2))
ggsave("figures/suppl_phylo_varcomp.pdf")

## Variance partitioning: niche vs migclass/body size

## Check spp overlap with mig distances from La Sorte
common_names <- clim %>%
  left_join(hab, by = c("species_code" = "spp")) %>%
  left_join(dplyr::select(diet, aou, shannonE_diet), by = c("aou")) %>%
  left_join(range) %>%
  left_join(select(poptrend, AOU, Trend), by = c("aou" = "AOU")) %>%
  filter(!is.na(mean_area) & !is.na(ssi)) %>%
  left_join(ro_correlates, by = c("aou" = "AOU")) %>%
  left_join(range_overlap) %>%
  select(english_common_name, aou, species_code, migclass)

mig_files <- list.files("/Users/gracedicecco/git/photoperiod-master/data/eBird")
mig_names <- data.frame(file = mig_files) %>%
  mutate(spp_common = word(file, 1, 1, sep =fixed('.'))) %>%
  mutate(english_common_name = gsub("_", " ", spp_common)) 
                        
join_mig <- common_names %>%
  mutate(english_common_name = gsub("'", " ", english_common_name)) %>%
  left_join(mig_names)%>%
  mutate(distance = ifelse(is.na(file), 0, 1))
table(join_mig$migclass, join_mig$distance)
# neotrop: 17 no, 31 yes
# short: 41 no, 1 yes
# resid: 16 no, 0 yes

varpart_vars <- all_vars %>%
  left_join(mig_unnest) %>%
  na.omit()

varpart_trend <- varpart(varpart_vars$Trend, ~ shannonE_diet, 
                       ~ mig_dist_m, 
                       ~ logMass,
                       data = varpart_vars)

varpart_occ <- varpart(varpart_vars$mean_occ, ~ ssi, 
                       ~ mig_dist_m, 
                       ~ logMass,
                       data = varpart_vars)

varpart_area <- varpart(varpart_vars$mean_area, ~ climate_vol, 
                        ~ mig_dist_m, 
                        ~ logMass,
                        data = varpart_vars)
