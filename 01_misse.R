setwd("~/Desktop/rhamnaceae_ploidy/")
rm(list=ls())

source("00_utility_functions.R")

phy <- read.tree("01-Pomaderris-example-spatial/RHAM2.215t.TreePL.LABELED.tree")
plot(phy, show.tip.label = F)
f = 0.86 #?
phy <- multi2di(phy)

max.param = max(4, round(ape::Ntip(phy)/10))
possible.combos = generateMiSSEGreedyCombinations(max.param=max.param, vary.both=TRUE, fixed.eps.tries=NA)

save.file = "Pomaderris_run.Rsave"
stop.deltaAICc = 2
n.cores = 4
chunk.size = 5

model.set = MiSSEGreedy(phy=phy, # the phylogeny 
                        f=f, # sampling fraction
                        possible.combos=possible.combos, # possible combinations of models
                        save.file=save.file, # the name of the file to save
                        stop.deltaAICc=stop.deltaAICc, # the deltaAIC to stop running
                        n.cores=n.cores, # number of cores
                        chunk.size=chunk.size, # number of models to run at each time
                        sann=FALSE) # IMPORTANT IMPORTANT IMPORTANT - This argument is here set to F for speed, 
# but it should be set to TRUE for your actual runs. See above. 

model.set_pruned <- PruneRedundantModels(model.set)

model.recons <- as.list(1:length(model.set_pruned))
for (model_index in 1:length(model.set_pruned)) {
  nturnover <- length(unique(model.set_pruned[[model_index]]$turnover))
  neps <- length(unique(model.set_pruned[[model_index]]$eps))
  model.recons[[model_index]] <- hisse::MarginReconMiSSE(phy = model.set_pruned[[model_index]]$phy, f = 0.86, hidden.states = max(c(nturnover, neps)), 
                                                         pars = model.set_pruned[[model_index]]$solution, fixed.eps=model.set_pruned$fixed.eps , 
                                                         AIC = model.set_pruned[[model_index]]$AIC, root.type = "madfitz",n.cores=n.cores)   
}


tip.rates <- GetModelAveRates(model.recons, type = "tips")
write.csv(tip.rates, file="misse_results/Pomaderris_tip_rates.csv", row.names=F)

#########################################################################
# Recommended to save all files:
save(model.set_pruned, 
     model.recons, 
     tip.rates, 
     possible.combos,
     file="misse_results/Pomaderris_results.Rsave")

#load("misse_results/Pomaderris_results.Rsave")
#########################################################################
# You can also visualize the rates evolving in the tree with the following command, 
# though rates through time *should not* be interpreted literally. The painting
# is just to get a "feeling" for the model.

pdf("misse_results/Pomaderris_misse.pdf", width=6, height=15)
rates <- c("net.div","speciation","turnover","extinction","extinct.fraction")
for(rate_index in 1:length(rates)){
  painted.tree <- hisse::plot.misse.states(x = model.recons, 
                                           rate.param = rates[rate_index], type = "phylo", show.tip.label = T, fsize=0.5)   
  title(main=rates[rate_index])
}

dev.off()
