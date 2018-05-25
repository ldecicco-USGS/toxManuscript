library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

source(file = "data_setup.R")

chem_site <- tox_list[["chem_site"]]

dir.create(file.path("plots"), showWarnings = FALSE)
pdf("plots/class_stacks.pdf", width = 11, height = 9)
i <- 1
for(class in unique(chemicalSummary$Class)){
  grid.newpage()
  sub_class <- filter(chemicalSummary, Class %in% class)
  
  cleaned_class <- tolower(class)
  
  if(cleaned_class == "pahs"){
    cleaned_class <- "PAHs"
  }
  
  fancyTitle <- paste0("Figure SI-4",LETTERS[i],
": Maximum exposure-activity ratio values by site for chemical class
", cleaned_class," in Great Lakes tributaries, 2010-2013.")
  i <- i+1
  upperPlot <- plot_tox_stacks(sub_class, 
                               chem_site, 
                               title = fancyTitle,
                               category = "Chemical")

  gb <- ggplot2::ggplot_build(upperPlot)
  gt <- ggplot2::ggplot_gtable(gb)
  
  gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
  grid::grid.draw(gt)
}
dev.off()
