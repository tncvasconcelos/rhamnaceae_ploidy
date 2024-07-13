
Description of folders: 
 
- **data/** 

Time callibrated phylogenetic tree and filtered distribution points for species of Pomaderreae included in the study. Also includes a template map used for mapping the distribution of the group in 

- **distribution_results/** 

Species distribution models results. Also includes save files for species-richness and mean net-diversification maps.

- **maps/** 

Species-richness and mean net-diversification rates maps following the distribution of Pomaderreae (as inferred using species distribution modeling).

- **misse_results/**

Files results from the MiSSEGreedy analysis. Includes list of models and their AICc, reconstruction files, tip-rates table and reconstruction plots of different diversification rate parameters along the tree (for visualization purposes only).


----
Description of scripts:

> 00_utility_functions.R

Functions used to model and map species distributions. Most functions in this script are internal of the GetRanges function.

> 01_misse.R

Pipeline to run MiSSEGreedy.

> 02_dist_models.R

Running GetRanges function.

> 03_species_richness.R

Calculates species-richness for Pomaderreae and Pomaderreae minus Pomaderris by overlapping binary layers and calculating their sum.

> 04_tip_rate_spatial.R

Calculates mean local net-diversification rates through the distribution of Pomaderreae by overlapping binary layers and assigning tip-rates to the distribution of each species.

> 05_map_plots.R

Plots maps for species-richness, mean net-diversification rates, and residuals of a linear regression between species-richness and net-diversification.
