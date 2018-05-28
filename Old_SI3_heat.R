library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

# Determine which chemicals are highest priority
EARthresh <- 0.001
siteThresh <- 10

SiteChemsCount <- filter(chemicalSummary,EAR>0) %>%
  group_by(shortName,date,chnm) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(shortName,chnm) %>%
  summarize(EARsummax = max(EARsum)) %>%
  filter(EARsummax > EARthresh) %>%
  group_by(chnm) %>%
  summarize(numSites = n_distinct(shortName)) %>%
  filter(numSites >= siteThresh)

priorityChems <- as.character(SiteChemsCount$chnm)
priorityChemRows <- which(chemicalSummary$chnm %in% priorityChems)


# NEED TO HIGHLIGHT THE PRIORITY CHEMICALS

heat <- plot_tox_heatmap(chemicalSummary,
                 tox_list$chem_site,
                 category = "Chemical")

dir.create(file.path("plots"), showWarnings = FALSE)

ggsave(heat, filename = "plots/SI3_Chem_heat.png", width = 10, height = 6)

#SI 5 no longer needed
# bio_heat <- plot_tox_heatmap(chemicalSummary,
#                  tox_list$chem_site,
#                  category = "Biological",
#                  manual_remove = "Undefined")
# 
# ggsave(bio_heat, filename = "plots/SI5_Bio_heat.png", width = 10, height = 6)
