rm(list=ls())

# Fig SI 2:
source("SI2_detection_levels.R", print.eval=TRUE)
rm(list=ls())

# Figs SI 3 & 6:
source("si_figs.R", print.eval=TRUE)
rm(list=ls())


# Fig SI 5:
source("SI5_chem_counts.R")
rm(list=ls())

# SI 4:
source("SI4_class_stacks.R")
rm(list=ls())


filenms <- c("SI2_detection_levels.png", "si3.png","SI5_chem_counts.png","si6.png")
file_path <- "./plots"

source("merge.png.files.R")
merge.png.pdf(pdfFile = "plots/SI_figs.pdf",pngFiles = filenms,file_path,deletePngFiles = FALSE)

SI_fig_captions <- character()
SI_fig_captions[1] <- "Figure SI-1. Location of sampling locations, watershed boundaries, and watershed land-uses (originally published in Baldwin et al., 2016). Map IDs are defined in Table SI-1."

SI_fig_captions[2] <- "Figure SI-2. Exposure activity ratios at the minimum detection level of each individual chemical measured in water samples at individual Great Lakes tributaries monitored from 2010-2013 for associated ToxCast endpoints."

SI_fig_captions[3] <- "Figure  SI-3. Maximum exposure activity ratios (EARs) for individual chemicals, summations within chemical classes, and summations for all chemicals measured in water samples at individual Great Lakes tributaries monitored from 2010-2013. "

SI_fig_captions[4] <- "Figure SI−4A-M: Maximum exposure−activity ratio values by site for chemical classes in Great Lakes tributaries, 2010−2013."

SI_fig_captions[5] <- "Figure SI-5. Number of sites with at least one sample that resulted in an exposure activity ratio > 10-3 for individual chemicals measured in water samples at individual Great Lakes tributaries monitored from 2010-2013."

SI_fig_captions[6] <- "Figure  SI-6. Maximum exposure activity ratios (EARs) by site for individual Adverse Outcome Pathways (AOPs), summations within chemical classes, and summations for all chemicals measured in water samples at individual Great Lakes tributaries monitored from 2010-2013."

SI_fig_captions[7] <- "Figure SI-7 Maximum exposure activity ratios (EARs) per site by Adverse Outcome Pathway (AOP) for mixtures of 2-, 3-, and 4-chemicals that occurred concurrently in samples with AOP EAR summations greater than 10-3 from four or more sites during the study period."

write(SI_fig_captions,file="plots/SI_fig_captions.txt",sep="")
