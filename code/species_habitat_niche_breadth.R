### Species habitat niche breadth

### Libraries

library(tidyverse)
library(purrr)
library(tmap)
library(sf)

### Plotting theme

theme_set(theme_classic())

### Data directories
bigdata <- "bigdata/"
rawdata <- "raw_data/"
derived_data <- "derived_data/"

## BioArk directory
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

## species range maps directory
range_dir <- paste0(bioark, "/hurlbertlab/GIS/birds/All/All/")

### Species list 

species_list <- read.csv(paste0(rawdata, "species_list.csv"), stringsAsFactors = F)

# Match BBS taxonomy with breeding range polygon taxonomy
fix_mismatch <- read.csv(paste0(derived_data, "fix_breedingrange_genus_mismatch.csv"), stringsAsFactors = F) %>%
  dplyr::select(-file, -spp_name) %>%
  filter(!is.na(old_genus)) %>%
  mutate(new_binomial = paste(old_genus, species))

bbs_spp <- species_list %>%
  mutate(binomial = paste(genus, species, sep = " ")) %>%
  filter(!grepl("unid.", english_common_name), !grepl("hybrid", english_common_name)) %>%
  mutate_at(c("binomial"), ~case_when(grepl("Colaptes auratus", .) ~ "Colaptes auratus",
                                      grepl("Junco hyemalis", .) ~ "Junco hyemalis",
                                      grepl("Setophaga coronata", .) ~ "Dendroica coronata",
                                      TRUE ~ .)) %>%
  left_join(fix_mismatch) %>%
  mutate(matched_name = ifelse(!is.na(old_genus), new_binomial, binomial),
         matched_filename = gsub(" ", "_", matched_name)) %>%
  filter(!is.na(matched_name))

range_files <- data.frame(file = list.files(range_dir)) %>%
  filter(grepl(".shp", file)) %>%
  mutate(spp_name = word(file, 1, 2, sep = "_"),
         file_binomial = gsub("_", " ", spp_name)) 

# Match with eBird codes

spp_list <- range_files %>%
  right_join(bbs_spp, by = c("file_binomial" = "matched_name")) %>%
  left_join(ebirdst::ebirdst_runs, by = c("english_common_name" = "common_name"))

### Function: extract centroid from breeding range polygons

range_centroid <- function(file) {
  br <- read_sf(paste0(range_dir, file))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
    st_transform(4326)
  
  cent <- st_centroid(breeding_range)
  
  coords <- st_coordinates(cent)
  
  return(coords)
}

### read in ebird occurrence averaged by land cover type
### processed on longleaf HPC

occ_lc <- data.frame(files = list.files(paste0(bigdata, "occ_lc"))) %>%
  group_by(files) %>%
  nest() %>%
  mutate(spp = word(files, 1, sep = "_")) %>%
  mutate(data = map(files, ~read_csv(paste0(bigdata, "occ_lc/", .)))) %>%
  unnest(cols = c("data")) %>%
  filter(!is.na(lc)) %>%
  group_by(spp) %>%
  mutate(nclass = n_distinct(lc))

## Calculate SSI for each species

ssi <- occ_lc %>%
  group_by(spp) %>%
  summarize(ssi = sd(meanOcc)/mean(meanOcc))

## Map of SSI vs range centroid for each species

# calculate species range centroids
cents <- spp_list %>%
  filter(!is.na(file)) %>%
  mutate(centroid = map(file, ~range_centroid(.)))

ssi_centroids <- cents %>%
  mutate(x = map_dbl(centroid, ~{
    c <- .
    df <- data.frame(c)
    
    mean(df$X, na.rm = T)
  }),
  y = map_dbl(centroid, ~{
    c <- .
    df <- data.frame(c)
    
    mean(df$Y, na.rm = T)
  })) %>%
  dplyr::select(-centroid) %>%
  left_join(ssi, by = c("species_code" = "spp")) %>%
  filter(!is.na(ssi)) %>%
  st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326)

na_map <- read_sf(paste0(rawdata, "ne_50m_admin_1_states_provinces_lakes.shp")) %>%
  filter(sr_adm0_a3 ==  "USA" | sr_adm0_a3 == "CAN") %>%
  filter(iso_3166_2 != "US-HI")

ssi_map <- tm_shape(na_map) + tm_polygons() +
  tm_shape(ssi_centroids) + tm_dots(col = "ssi", palette = "YlGn", size = 0.4, title = "SSI - habitat")
tmap_save(ssi_map, "figures/ssi_habitat_niche_map.pdf")

#### Old code #####

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


