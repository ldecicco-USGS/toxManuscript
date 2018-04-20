library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

source(file = "data_setup.R")

chem_site <- tox_list[["chem_site"]]

pdf("class_stacks.pdf", width = 11, height = 9)
for(class in unique(chemicalSummary$Class)){
  grid.newpage()
  sub_class <- filter(chemicalSummary, Class %in% class)
  
  upperPlot <- plot_tox_stacks(sub_class, 
                               chem_site, 
                               title = class,
                               category = "Chemical")

  gb <- ggplot2::ggplot_build(upperPlot)
  gt <- ggplot2::ggplot_gtable(gb)
  
  gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
  grid::grid.draw(gt)
}
dev.off()
