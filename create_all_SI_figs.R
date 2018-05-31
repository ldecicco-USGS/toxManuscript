rm(list=ls())

# Fig SI 2:
source("SI2_detection_levels.R", print.eval=TRUE)
rm(list=ls())

# Figs SI 3 & 6:
source("si_figs.R", print.eval=TRUE)
rm(list=ls())


# Fig SI 2:
source("SI5_chem_counts.R")
rm(list=ls())

# SI 4:
source("SI4_class_stacks.R")
rm(list=ls())


filenms <- c("SI2_detection_levels.png", "si3.png","SI5_chem_counts.png","si6.png")
file_path <- "./plots"

source("merge.png.files.R")
merge.png.pdf(pdfFile = "plots/SI_figs.pdf",pngFiles = filenms,file_path,deletePngFiles = FALSE)
