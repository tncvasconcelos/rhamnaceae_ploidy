# rm(list=ls())
# setwd("~/Desktop/rhamnaceae_ploidy/")

source("00_utility_functions.R")

all_dist_points <- read.csv("data/C-RHAM_GEO_ALL_2023-03-02-extra-col.csv")

all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Longitude))
all_dist_points <- subset(all_dist_points, !is.na(all_dist_points$Latitude))

test_run <- GetRanges(points=all_dist_points, species="Species", lat="Latitude", lon="Longitude", threshold="", buffer=25, res=10)
saveRDS(test_run, file="distribution_models_10th_percentile_thrhld.Rdata")

