### Species habitat niche breadth

### Libraries

library(tidyverse)

### Plotting theme

theme_set(theme_classic())

### Source function to get eBirdst PIs for a species

source("code/get_ebird_pis.R")

### Species list

species <- read_csv("raw_data/species_list.csv")

ebird_bbs <- ebirdst::ebirdst_runs %>% 
  left_join(species, by = c("common_name" ="english_common_name")) %>%
  filter(!is.na(french_common_name))

### For each species, create a PI csv and save to /pis
# easpho & easblu bug - skip

p <- "/bigdata/pis/"

for(s in spp_left$species_code) {
  
  get_breedingrange_pis(path = p, spp_code = s)
  
  print(paste0(Sys.time(), " species complete: ", s))
  
}

## Calculate Shannon Index of model land cover PI for each species

these_predictors <- ebirdst::ebirdst_predictors

lc_pred <- these_predictors %>%
  filter(!is.na(lc_class))

spp_pi <- data.frame(spp_code = word(list.files("bigdata/pis"), 1, 1, sep = "_"),
                    filepath = list.files("bigdata/pis"))

spp_files <- spp_pi %>%
  group_by(spp_code, filepath) %>%
  nest() %>%
  mutate(data = map(filepath, ~read_csv(paste0(getwd(),"/bigdata/pis/",.))))

spp_hab_index <- spp_files %>%
  mutate(med_pi = map(data, ~{
    df <- .
    
    # remove spurious large values per documentation
    
    df_subs <- df %>%
      pivot_longer(names_to = "Predictor", values_to = "PI", cols = everything()) %>%
      filter(PI < quantile(PI, probs = 0.98, na.rm = T))
    
    df_subs %>%
      group_by(Predictor) %>%
      summarize(PI = median(PI, na.rm = T))
    
    }),
    lc_evenness = map_dbl(med_pi, ~{
      df <- .
      
      lc_only <- df %>%
        filter(Predictor %in% lc_pred$lc_class_label) %>%
        na.omit()
      
      hmax <- log(nrow(lc_only))
      totalpi <- sum(lc_only$PI)
      
      no_zeroes <- lc_only %>% 
        filter(PI > 0)
      
      h <- -sum( (no_zeroes$PI/totalpi) * log(no_zeroes$PI/totalpi))
      
      h/hmax
      
    }),
    lc_dominance = map_dbl(med_pi, ~{
      df <- .
      
      lc_only <- df %>%
        filter(Predictor %in% lc_pred$lc_class_label) %>%
        na.omit()
      
      totalpi <- sum(lc_only$PI)
      
      maxpi <- max(lc_only$PI)
      
      maxpi/totalpi
    }),
    lc_dominance_class = map_chr(med_pi, ~{
      df <- .
      
      lc_only <- df %>%
        filter(Predictor %in% lc_pred$lc_class_label) %>%
        na.omit()
      
      maxpi <- max(lc_only$PI)
      
      lc_only$Predictor[lc_only$PI == maxpi]
    })) 

spp_hab_breadth <- spp_hab_index %>%
  select(-data, -med_pi) %>%
  left_join(ebird_bbs, by = c("spp_code" = "species_code")) %>%
  ungroup() %>%
  select(spp_code, common_name, lc_evenness, lc_dominance, lc_dominance_class)
write.csv(spp_hab_breadth, "derived_data/spp_habitat_niche_ebird.csv", row.names = F)

# Plot histograms of Shannon E, dominance metrics

ggplot(spp_hab_index, aes(x = lc_evenness)) + geom_histogram(col = "white") +
  labs(x = "Shannon Evenness - Habitat", y = "Species")
ggsave("figures/habitat_niche_shannonE_hist.pdf")

ggplot(spp_hab_index, aes(x = lc_dominance)) + geom_histogram(col = "white") + 
  labs(x = "Dominance - Habitat", y = "Species")
ggsave("figures/habitat_niche_dominance_hist.pdf")

ggplot(spp_hab_index, aes(x = lc_dominance_class)) + geom_histogram(stat = "count", col = "white") +
  coord_flip() + labs(x = "Dominant class - Habitat", y = "Species")
ggsave("figures/habitat_niche_dominance_class_counts.pdf")

# Plot density of median PI values for each LC predictor class

predictor_importance <- spp_hab_index %>%
  select(-data) %>%
  unnest(cols = c("med_pi")) %>%
  filter(Predictor %in% lc_pred$lc_class_label)

plot_lc_density <- function(data, Predictor) {
  plot <- ggplot(data, aes(x = PI)) + 
    geom_density(fill = "gray") + labs(x = "PI", title = Predictor)
  return(plot)
}

plot_list <- predictor_importance %>%
  group_by(Predictor) %>%
  nest() %>%
  pmap(plot_lc_density)

multi_page <- gridExtra::marrangeGrob(plot_list, nrow = 2, ncol = 3)
ggsave("figures/habitat_niche_landcover_density.pdf", multi_page, units = "in", height = 6, width = 10)

## Plot PIs and Shannon Index for each species in multi-page PDF

plot_hab <- function(common_name, lc_evenness, lc_dominance, med_pi, ...) {
  pi_long <- med_pi %>%
    filter(Predictor %in% lc_pred$lc_class_label) 
  
  theme_set(theme_classic(base_size = 15))
  
  ggplot(pi_long, aes(x = fct_reorder(Predictor, PI), y = PI)) + geom_point() +
    labs(title = paste0(common_name),
         x = NULL, y = "Relative PI") +
    annotate(geom = "text", x = 3, y = 0.5*max(pi_long$PI), label = paste0("E = ", round(lc_evenness, 2)), size = 6) +
    annotate(geom = "text", x = 1, y = 0.5*max(pi_long$PI), label = paste0("D = ", round(lc_dominance, 2)), size = 6) +          
    coord_flip()
}

spp_plots <- spp_hab_index %>%
  left_join(ebird_bbs, by = c("spp_code" = "species_code")) %>%
  arrange(aou) %>%
  pmap(plot_hab)

spp_multi <- gridExtra::marrangeGrob(spp_plots, nrow = 1, ncol = 1)
ggsave("figures/habitat_niche_spp_landcover_evenness.pdf", spp_multi, units = "in", height = 6, width = 8)


