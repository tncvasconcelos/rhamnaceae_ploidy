setwd("~/Desktop/rhamnaceae_ploidy/")
#rm(list=ls())
library(maptools)
data("wrld_simpl")

# Load models again
test_run <- readRDS("distribution_models.Rdata")

# Removing models that didn't work
to_remove<-c()
for(i in 1:length(test_run)) {
  if(length(test_run[[i]])==1){
    to_remove <- c(to_remove, i)
  }
}
test_run[to_remove] <- NULL

# Load tip rates
tip_rates <- read.csv("misse_results/Pomaderris_tip_rates.csv")

library(raster)
# Organize plot of species richness
ranges <- unlist(lapply(test_run, "[[", "range"))
template.map <- readRDS("data/template.map.Rdata")
# Cut for Australia to make things faster
focal_extent <- c(100,180,-50,0)
template.map <- crop(template.map, focal_extent)
#ranges <- lapply(ranges, crop, focal_extent)

names(ranges) <- gsub(" ","_",names(ranges))
ranges <- subset(ranges, names(ranges) %in% tip_rates$taxon)

tmp.raster.list <- list()
for (i in 1:length(ranges)) {
  one_label <- names(ranges)[i]
  r1 <- ranges[[i]]
  r1 <- raster::resample(r1, template.map)
  r1[is.na(r1)] <- 0
  r1[which(r1[]==1)] <- tip_rates$net.div[tip_rates$taxon==one_label]
  tmp.raster.list[[i]] <- raster::mask(r1, template.map)
  print(i)
}
netdiv_map <- raster::calc(raster::stack(tmp.raster.list), mean)

saveRDS(netdiv_map, file="netdiv_map.Rdata")

pdf("netdiv_map.pdf")
raster::plot(netdiv_map)
dev.off()

######
# netdiv without Pomaderris

ranges_woPomaderris <- subset(ranges ,!grepl("Pomaderris", names(ranges)))
tmp.raster.list <- list()
buffer_mask=10
for (i in 1:length(ranges_woPomaderris)) {
  one_label <- names(ranges_woPomaderris)[i]
  r1 <- ranges_woPomaderris[[i]]
  r1 <- rasterToPoints(r1)
  if(nrow(r1)>0) {
    circle_around_point <- dismo::circles(r1, d=buffer_mask*1000, lonlat=TRUE)
    circle_around_point <- polygons(circle_around_point)
    r1 <- rasterize(circle_around_point, template.map)
    r1 <- raster::resample(r1, template.map)
    r1[is.na(r1)] <- 0
    r1[which(r1[]==1)] <- tip_rates$net.div[tip_rates$taxon==one_label]
    r0 <- raster::mask(r1, template.map)
    r0[is.na(template.map)] <- NA
    tmp.raster.list[[i]] <- r0
    print(i)    
  }
}
tmp.raster.list[which(unlist(lapply(tmp.raster.list, is.null)))] <- NULL
netdiv_map <- raster::calc(raster::stack(tmp.raster.list), mean)
saveRDS(netdiv_map, file="netdiv_map_woPomaderris.Rdata")

pdf("sp_rich_maps/netdiv_map_woPomaderris2.pdf")
raster::plot(netdiv_map)
dev.off()
