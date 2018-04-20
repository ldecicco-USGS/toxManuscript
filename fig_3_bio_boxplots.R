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

grid::grid.draw(gt)
