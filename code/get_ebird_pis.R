### Function: get eBird PIs for a species during breeding season, all of North America

library(ebirdst)

# Inputs:
# path = where to save PIs csv
# spp_code = ebird species code (6 letters)
get_breedingrange_pis <- function(path, spp_code) {
  sp_path <- ebirdst_download(species = spp_code, tifs_only = F)
  
  pis <- load_pis(sp_path)
  
  # All of North America, only June (breeding)
  lp_extent <- ebirdst_extent(c(xmin = -180, xmax = 50, ymin = 25, ymax = 70),
                              t = c("2018-06-01", "2018-06-30"))
  
  ### From source: getting underlying df from plot_pis fn
  
  these_predictors <- ebirdst::ebirdst_predictors
  
  # subset
  pis_ext <- ebirdst_subset(pis, ext = lp_extent)
  pis_ext <- pis[, these_predictors$predictor_tidy]
  
  # if aggregating by cover class aggregate the fragstats metrics
  # find landcover classes
  lc <- convert_classes(names(pis_ext), by_cover_class = TRUE,
                        pretty = T)
  lc_groups <- unique(lc)
  
  # aggregate over classes
  m <- matrix(nrow = nrow(pis_ext), ncol = length(lc_groups))
  colnames(m) <- lc_groups
  
  for (i in lc_groups) {
    if (sum(lc == i) == 1) {
      m[, i] <- pis_ext[, lc == i]
    } else {
      m[, i] <- apply(pis_ext[, lc == i], 1, FUN = mean, na.rm = TRUE)
    }
  }
  
  pis_lc <- as.data.frame(m, stringsAsFactors = FALSE)
  
  # write PIs data frame
  write.csv(pis_lc, paste0(getwd(), path, spp_code, "_breeding_pis.csv"), row.names = F)
  
  # remove species file to save space (~1.5 gb per species)
  unlink(sp_path, recursive = T)
}