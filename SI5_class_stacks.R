library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

source(file = "data_setup.R")
source(file = "plot_tox_stacks_manuscript.R")

chem_site <- tox_list[["chem_site"]]

dir.create(file.path("plots"), showWarnings = FALSE)

i <- 1
for(class in unique(chemicalSummary$Class)){

  sub_class <- filter(chemicalSummary, Class %in% class)
  
  cleaned_class <- tolower(class)
  
  if(cleaned_class == "pahs"){
    cleaned_class <- "PAHs"
  }
  
  y_label <- toxEval:::fancyLabels("Chemical", FALSE, TRUE, FALSE, sep = TRUE)
  y_label[["y_label"]] <- bquote("max" ~  group("(", sum(" " ~ EAR["[" * i * "]"]), ")")["[" * j * "]"]) 
  y_label[["caption"]] <- gsub(", k = sites","",y_label[["caption"]])
  
  fancyTitle <- bquote(atop(bold("Figure SI-5"~.(LETTERS[i])~":") ~
"Maximum exposure-activity ratio values by site for chemical class" ~ .(cleaned_class) ~ " in Great Lakes tributaries, 2010-2013. (" ~
  .(y_label[["caption"]]) ~ ")"))
  
  splot <- plot_tox_stacks_manuscript(sub_class, 
                               chem_site, 
                               category = "Chemical",
                               caption = fancyTitle)
  ggsave(splot, filename = paste0("plots/SI-5-",i,".pdf"), width = 11, height = 5)
  i <- i+1
}


