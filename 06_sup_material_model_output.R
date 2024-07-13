#setwd("~/Desktop/rhamnaceae_ploidy/")
#rm(list=ls())
library(maptools)
library(raster)
data("wrld_simpl")

# Load models again
test_run <- readRDS("distribution_models_10th_percentile_thrhld.Rdata")
# creating table with AUCs and thresholds to go into the sup material
result_sup_material <- as.data.frame(matrix(nrow=length(test_run), ncol=5))
for(i in 1:length(test_run)) {
  one_species <- test_run[[i]]
  result_sup_material[i,1] <- one_species$species_name
  result_sup_material[i,2] <- one_species$n_points
  if("threshold" %in% names(one_species)) {
    result_sup_material[i,3] <- paste(one_species$predictors, collapse = " ")
    result_sup_material[i,4] <- one_species$threshold
    result_sup_material[i,5] <- one_species$auc
  } else {
    result_sup_material[i,3] <- NA
    result_sup_material[i,4] <- NA
    result_sup_material[i,5] <- one_species$note
  }
}
colnames(result_sup_material) <- c("Species","n points SDM","predictors","10th percentile threshold","AUC")
write.csv(result_sup_material, "SDM_sup_material.csv", row.names = F)
