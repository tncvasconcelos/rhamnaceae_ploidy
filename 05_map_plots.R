#setwd("~/Desktop/rhamnaceae_ploidy/")
#rm(list=ls())
library(RColorBrewer)
library(maptools)
library(raster)
data("wrld_simpl")
source("00_utility_functions.R")

netdiv_map <- readRDS("netdiv_map.Rdata")
sprichness_map <- readRDS("species_richness_map.Rdata")
res_test <- GetResDistribution(raster1=sprichness_map, raster2=netdiv_map)

australia <- subset(wrld_simpl, NAME=="Australia")

pdf("maps.pdf", height=5, width=10)
# plots
par(mfrow=c(1,3))
par(oma=c(1,1,1,1))
#par(mar = c(2, 5, 2, 5))
# sp_rich
pal <- hcl.colors(30, palette = "RdYlBu", alpha = 1)
sprichness_map_aus <- mask(sprichness_map, australia)
raster::plot(sprichness_map_aus, col=rev(pal), xlim=c(110,160))
title(main="species richness")
plot(australia, add=T)

# net_div
netdiv_map_aus <- mask(netdiv_map, australia)
raster::plot(netdiv_map_aus, col=rev(pal), xlim=c(110,160))
title(main="net diversification")
plot(australia, add=T)

# regression
neg_values <- seq(from=min(res_test[], na.rm=T), to=0, length.out=9)
pos_values <- seq(from=0, to= max(res_test[], na.rm=T), length.out=9)
breakpoints <- unique(round(c(neg_values,0, pos_values), 2))
colors <- c(rev(RColorBrewer::brewer.pal(7, "Reds")),"#FFFFFF", RColorBrewer::brewer.pal(8, "Blues"))
plot(res_test, breaks=breakpoints, col = colors, xlim=c(110,160))
title(main="species richness~net diversification (residuals)")
plot(australia, add=T)

dev.off()

# sp rich without Pomaderris
species_richness_map_woPomaderris <- readRDS("species_richness_map_woPomaderris.Rdata")
plot(species_richness_map_woPomaderris)

pdf("sprich_wo_Pomaderris.pdf")
pal <- hcl.colors(30, palette = "RdYlBu", alpha = 1)
species_richness_map_woPomaderris <- mask(species_richness_map_woPomaderris, australia)
raster::plot(species_richness_map_woPomaderris, col=rev(pal), xlim=c(110,160))
title(main="species richness wo Pomaderris")
plot(australia, add=T)

dev.off()

