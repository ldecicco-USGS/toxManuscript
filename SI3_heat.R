library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

heat <- plot_tox_heatmap(chemicalSummary,
                 tox_list$chem_site,
                 category = "Chemical")
