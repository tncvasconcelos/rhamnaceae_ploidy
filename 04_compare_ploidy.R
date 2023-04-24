setwd("~/Desktop/rhamnaceae_ploidy/")
# rm(list=ls())

tip_rates <- read.csv("misse_results/Pomaderris_tip_rates.csv")
ploidy <- read.csv("00-RHAM-Phylo-217t-PLOIDY-0203/RHAM_PLOIDY-match.Treepl.217t.0203.csv")
colnames(tip_rates)[which(colnames(tip_rates)=="taxon")] <- "species"
merged <- merge(tip_rates, ploidy, by="species")

pdf("ploidy_netdiv_misse.pdf")
boxplot(merged$net.div~merged$ploidy, ylab="net.div",xlab="ploidy")
dev.off()

