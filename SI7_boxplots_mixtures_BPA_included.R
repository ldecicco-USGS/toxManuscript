library(toxEval)
library(dplyr)
library(tidyr)
library(data.table)
#########################################################################################

source(file = "Table_SI7_fig4_Data_summaries_for_manuscript.R")

#This script requires running Table_SI7_fig4_Data_summaries_for_manuscript.R and saving the
#data frame "Num_sites_by_mixture.csv" before running these plots.

EAR_thresh <- 0.001
# ep_percent_thres <- 0.5

chemSummData_max <- chemicalSummary %>%
  filter(EAR > 0) %>%
  left_join(AOP_relevance, by="endPoint") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE)) 


EAR_thresh <- 0.00001


plot_dimensions <- list(c(0,0),c(2,2),c(3,2),c(1,2))
margins <- c(4,0.5,1,0)
outer_margins <- c(7,5,2,1)
axis_text_cex <- 0.6
title_text_cex <- 0.7

y_label <- bquote(EAR[SiteAOP])
ToxCast_ACC$chnm[!is.na(ToxCast_ACC$chnm) & ToxCast_ACC$chnm == "TDCPP"] <- "Tris(1,3-dichloro-2-propyl)phosphate"
###################################
i <- 2
filenm <- "plots/BPA_mixtures/SI7_mixtureBoxplots_2_chems_BPA.pdf"
pdf(filenm)

sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)

par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)

k <- 0
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites$chemVector[j],"\\|")[[1]]
  if("80-05-7" %in% CASnums) {
    k <- k+1
  nMixSites <- sub_Num_sites$numSites[j]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(k %in%  c(1,3),"s","n")
  
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  title(paste(chnms, collapse = "\n"), line=-1.5, cex.main = title_text_cex,adj=0.03)
  mtext(paste(nMixSites,"Sites"),side=3,line=0,cex=0.7)
  
  }
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
# mtext(paste(i,"-Compound Mixtures"),outer=TRUE)
mtext(bquote(.(y_label)),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label)  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                    "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013.")),outer=TRUE)

dev.off()

###################################
i <- 3
filenm <- "plots/BPA_mixtures/SI7_mixtureBoxplots_3_chems_BPA.pdf"
pdf(filenm)
sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)
par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)
k <- 0
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites$chemVector[j],"\\|")[[1]]
  if("80-05-7" %in% CASnums) {
    k <- k + 1
    nMixSites <- sub_Num_sites$numSites[j]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(k %in%  c(1,3,5),"s","n")
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  title(paste(chnms, collapse = "\n"), line=-2.5, cex.main = title_text_cex,adj=0.03)
  mtext(paste(nMixSites,"Sites"),side=3,line=0,cex=0.7)
  }
  
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
# mtext(paste(i,"-Compound Mixtures"),outer=TRUE)
mtext(bquote(.(y_label[["y_label"]])),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label)  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                         "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013.")),outer=TRUE)

dev.off()

################################
i <- 4
filenm <- "plots/BPA_mixtures/SI7_mixtureBoxplots_4_chems_BPA.pdf"
pdf(filenm, height = 4)
sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)
par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites$chemVector[j],"\\|")[[1]]
  if("80-05-7" %in% CASnums) {
    nMixSites <- sub_Num_sites$numSites[j]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(j %in% 1,"s","n")
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  title(paste(chnms, collapse = "\n"), line=-3, cex.main = title_text_cex,adj=0.03)
  mtext(paste(nMixSites,"Sites"),side=3,line=0,cex=0.7)
  
  
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
# mtext(paste(i,"-Compound Mixtures"),outer=TRUE)
mtext(bquote(.(y_label[["y_label"]])),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label)  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                         "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013.")),outer=TRUE)

}
dev.off()
