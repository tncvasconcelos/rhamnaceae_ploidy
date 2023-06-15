
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

# species richness without Pomaderris

ranges_woPomaderris <- subset(ranges ,!grepl("Pomaderris", names(ranges)))
tmp.raster.list <- list()
for (i in 1:length(ranges_woPomaderris)) {
  r1 <- ranges_woPomaderris[[i]]
  r1 <- raster::resample(r1, template.map)
  r1[is.na(r1)] <- 0
  tmp.raster.list[[i]] <- raster::mask(r1, template.map)
  print(i)
}
names(tmp.raster.list) <- names(ranges_woPomaderris)
sprichness_map <- raster::calc(raster::stack(tmp.raster.list), sum)
saveRDS(sprichness_map, file="species_richness_map_woPomaderris.Rdata")

pdf("species_richness_map_woPomaderris.pdf")
raster::plot(sprichness_map)
dev.off()