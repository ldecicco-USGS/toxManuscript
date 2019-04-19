# Fig SI 2:
source("SI2_detection_levels.R", print.eval=TRUE)
rm(list=ls())

# Fig SI 3:
source("SI3_chem_counts.R", print.eval=TRUE)
rm(list=ls())

# Fig SI 4:
source("SI4_heat_map.R", print.eval=TRUE)
rm(list=ls())

# Fig SI 5:
source("SI5_class_stacks.R", print.eval=TRUE)
rm(list=ls())

# Fig SI 6:
source("SI6_AOP_heat.R", print.eval=TRUE)
rm(list=ls())

# Fig SI 7:
source("SI7_boxplots_mixtures_For_AOP_networks_fig.R", print.eval=TRUE) #Changed this to include only those in fig 6
rm(list=ls())

# library(staplr)
# 
# working_dir <- getwd()
# output_file <- "SI_figs.pdf"
# plots_dir <- file.path(working_dir,"plots")
# #setwd(plots_dir)
# staple_pdf(input_directory = plots_dir, input_files = NULL, output_filepath = output_file)
# #setwd(working_dir)
