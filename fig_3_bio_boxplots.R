library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

bioPlot <- plot_tox_boxplots(chemicalSummary, 
                             category = "Biological", 
                             manual_remove = c("Transferase","Undefined"))

gb <- ggplot2::ggplot_build(bioPlot)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel"] <- "off"

dir.create(file.path("plots"), showWarnings = FALSE)
png("plots/bio_boxplots.png", width = 1000, height = 800, res = 142)
grid::grid.draw(gt)
dev.off()