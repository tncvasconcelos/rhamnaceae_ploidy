setwd("~/Desktop/rhamnaceae_ploidy/")
rm(list=ls())

source("00_utility_functions.R")

all_dist_points <- read.csv("02-RHAM-spatial-all/C-RHAM_GEO_ALL_2023-03-02-extra-col.csv")

all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Longitude))
all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Latitude))

test_run <- GetRanges(points=all_dist_points, species="Species", lat="Latitude", lon="Longitude", threshold=0.75, buffer=25, res=2.5)
#saveRDS(test_run, file="distribution_models.Rdata")

test_run <- readRDS("distribution_models.Rdata")

# Removing models that didn't work
to_remove<-c()
for(i in 1:length(test_run)) {
  if(length(test_run[[i]])==1){
    to_remove <- c(to_remove, i)
  }
}
test_run[to_remove] <- NULL

library(raster)
# Organize plot of species richness
ranges <- unlist(lapply(test_run, "[[", "range"))
template.map <- readRDS("data/template.map.Rdata")
# Cut for Australia to make things faster
focal_extent <- c(100,180,-50,0)
template.map <- crop(template.map, focal_extent)
#ranges <- lapply(ranges, crop, focal_extent)

tmp.raster.list <- list()
for (i in 1:length(ranges)) {
  r1 <- ranges[[i]]
  r1 <- raster::resample(r1, template.map)
  r1[is.na(r1)] <- 0
  tmp.raster.list[[i]] <- raster::mask(r1, template.map)
  print(i)
}
names(tmp.raster.list) <- names(ranges)
sprichness_map <- raster::calc(raster::stack(tmp.raster.list), sum)
saveRDS(sprichness_map, file="species_richness_map.Rdata")

pdf("species_richness_map.pdf")
raster::plot(sprichness_map, file="species_richness_map.Rdata")
dev.off()
 

