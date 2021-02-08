### Species habitat niche breadth

### Libraries

library(tidyverse)

### Source function to get eBirdst PIs for a species

source("get_ebird_pis.R")

### Species list

species <- read_csv("species_list.csv")

ebird_bbs <- ebirdst_runs %>% 
  left_join(species, by = c("common_name" ="english_common_name")) %>%
  filter(!is.na(french_common_name))

### For each species, create a PI csv and save to /pis
# easpho & easblu bug - skip

p <- "/pis/"

for(s in spp_left$species_code) {
  
  get_breedingrange_pis(path = p, spp_code = s)
  
  print(paste0(Sys.time(), " species complete: ", s))
  
}

## Calculate Shannon Index of model land cover PI for each species

lc_pred <- these_predictors %>%
  filter(!is.na(lc_class))

spp_pi <- data.frame(spp_code = word(list.files("pis"), 1, 1, sep = "_"),
                    filepath = list.files("pis"))

spp_files <- spp_pi %>%
  group_by(spp_code, filepath) %>%
  nest() %>%
  mutate(data = map(filepath, ~read_csv(paste0(getwd(),"/pis/",.))))

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
      
    }))

ggplot(spp_hab_index, aes(x = lc_evenness)) + geom_histogram(col = "white") +
  labs(x = "Shannon Evenness - Habitat", y = "Species")
ggsave("figures/shannonE_habitat_niche_hist.pdf")

## Plot PIs and Shannon Index for each species in multi-page PDF

plot_hab <- function(common_name, lc_evenness, data, ...) {
  pi_long <- data %>%
    pivot_longer(names_to = "Predictor", values_to = "PI", cols = everything()) %>%
    filter(PI < quantile(PI, probs = 0.98, na.rm = T)) %>%
    group_by(Predictor) %>%
    mutate(med = median(PI))
  
  theme_set(theme_classic(base_size = 15))
  
  p <- ggplot(pi_long, aes(x = fct_reorder(Predictor, med), y = PI)) + geom_boxplot() +
    labs(title = paste0(common_name, " E = ", round(lc_evenness, 3)),
         x = NULL, y = "Relative PI") +
    coord_flip()
  
  print(p)
}

for(i in 1:8) {
  start <- 45*i - 44
  
  if(i != 8) {
    end <- 45*i }
  else { end <- nrow(spp_hab_index) }
  
  pdf(paste0(getwd(), "/figures/spp_landcover_evenness_", start, "-", end, ".pdf"), height = 6, width = 8)
  
  subs <- spp_hab_index[start:end, ]
  
  subs %>%
    left_join(ebird_bbs, by = c("spp_code" = "species_code")) %>%
    arrange(aou) %>%
    pmap(plot_hab)
  
  dev.off()
}



