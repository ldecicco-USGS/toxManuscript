library(toxEval)
library(dplyr)
library(tidyr)

source(file = "data_setup.R")

heat_map <- plot_tox_heatmap(chemicalSummary, tox_list$chem_site, category = "Chemical")

ggsave(heat_map, filename = "plots/SI4_heat_map.png", height = 9, width = 11)

