# Create all the figures....

# Fig 1:
source("fig1_combo_graph.R")
rm(list=ls())
# Fig 2:
source("fig_2_site_counts.R")
rm(list=ls())
# Fig 3:
source("fig_3_endpoint_boxplots.R")
rm(list=ls())
# Fig 4:
source("fig_4_aop_priority.R", print.eval=TRUE)
rm(list=ls())


filenms <- c("fig1_no_clip.png", "fig2_site_count.png","fig3_endpoint_boxplots.png","Fig4_aop_cow.png")
file_path <- "./plots"

source("merge.png.files.R")
merge.png.pdf(pdfFile = "plots/figs.pdf",pngFiles = filenms,file_path,deletePngFiles = FALSE)


Fig_captions <- character()
Fig_captions[1] <- "Figure SI-1. Location of sampling locations, watershed boundaries, and watershed land-uses (originally published in Baldwin et al., 2016). Map IDs are defined in Table SI-1."

Fig_captions[2] <- ""

Fig_captions[3] <- ""

Fig_captions[4] <- ""

Fig_captions[5] <- ""
