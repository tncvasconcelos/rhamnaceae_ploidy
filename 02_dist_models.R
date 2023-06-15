setwd("~/Desktop/rhamnaceae_ploidy/")
#rm(list=ls())

source("00_utility_functions.R")

all_dist_points <- read.csv("02-RHAM-spatial-all/C-RHAM_GEO_ALL_2023-03-02-extra-col.csv")

all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Longitude))
all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Latitude))

test_run <- GetRanges(points=all_dist_points, species="Species", lat="Latitude", lon="Longitude", threshold=0.5, buffer=25, res=10)
#saveRDS(test_run, file="distribution_models.Rdata")

