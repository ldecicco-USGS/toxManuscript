library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

heat <- plot_tox_heatmap(chemicalSummary,
                 tox_list$chem_site,
                 category = "Chemical")

dir.create(file.path("plots"), showWarnings = FALSE)

ggsave(heat, filename = "plots/SI3_Chem_heat.png", width = 10, height = 6)

#SI 5
bio_heat <- plot_tox_heatmap(chemicalSummary,
                 tox_list$chem_site,
                 category = "Biological",
                 manual_remove = "Undefined")

ggsave(bio_heat, filename = "plots/SI5_Bio_heat.png", width = 10, height = 6)
