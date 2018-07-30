library(toxEval)
library(dplyr)
library(tidyr)
library(data.table)
#########################################################################################
source("data_setup.R")
source(file = "data_setup.R")
source(file = "MakeTitles.R")

#This script requires running fig4_Data_summaries_for_manuscript.R and saving the
#data frame "Num_sites_by_mixture.csv" before running these plots.

EAR_thresh <- 0.001
# ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)
relevance$Relevant <- MakeTitles(relevance$Relevant)

relevance <- relevance %>%
  select(ID=AOP,Relevant,Rationale)


AOP_relevance <- left_join(AOP,relevance,by="ID")


Num_sites_by_mixture <- read.csv(file="SI_table7 Num_sites_by_mixture_temp.csv",stringsAsFactors = FALSE)

chemSummData_max <- chemicalSummary %>%
  filter(EAR > 0) %>%
  left_join(AOP_relevance, by="endPoint") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE)) 


EAR_thresh <- 0.00001


plot_dimensions <- list(c(0,0),c(3,4),c(2,3),c(1,3))
margins <- c(4,0.5,1,0)
outer_margins <- c(7,5,2,1)
axis_text_cex <- 0.6
title_text_cex <- 0.7

y_label <- toxEval:::fancyLabels("Biological", FALSE, TRUE, FALSE, TRUE)
y_label$caption <- "i = chemicals in an AOP, j = samples, k = sites"
###################################
i <- 2
filenm <- "plots/SI7_mixtureBoxplots_A.pdf"
pdf(filenm)
sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)
par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites[j,"chemVector"],"; ")[[1]]
  nMixSites <- sub_Num_sites[j,"numSites"]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$casn %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(j %in%  (0:4*plot_dimensions[[i]][2] +1),"s","n")
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",main=paste(chnms),
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  mtext(paste(nMixSites,"Sites"),side=3,line=-1,cex=0.7)
  
  
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
mtext(bquote(.(y_label[["y_label"]])),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label[["y_label"]])  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                    "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013 (" * italic(.(y_label$caption)) * ").")),outer=TRUE)

dev.off()

###################################
i <- 3
filenm <- "plots/SI7_mixtureBoxplots_B.pdf"
pdf(filenm)
sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)
par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites[j,"chemVector"],"; ")[[1]]
  nMixSites <- sub_Num_sites[j,"numSites"]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$casn %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(j %in%  (0:4*plot_dimensions[[i]][2] +1),"s","n")
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",main=paste(chnms),
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  mtext(paste(nMixSites,"Sites"),side=3,line=-1,cex=0.7)
  
  
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
mtext(bquote(.(y_label[["y_label"]])),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label[["y_label"]])  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                         "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013 (" * italic(.(y_label$caption)) * ").")),outer=TRUE)

dev.off()

################################
i <- 4
filenm <- "plots/SI7_mixtureBoxplots_C.pdf"
pdf(filenm, height = 4)
sub_Num_sites <- Num_sites_by_mixture %>%
  filter(nChems == i,numSites>=4)
par(mfrow=plot_dimensions[[i]],mar=margins,oma=outer_margins)
for(j in 1:dim(sub_Num_sites)[1]){
  CASnums <- strsplit(sub_Num_sites[j,"chemVector"],"; ")[[1]]
  nMixSites <- sub_Num_sites[j,"numSites"]
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$casn %in% CASnums),"chnm"])
  
  subChemSummary <- chemSummData_max %>%
    filter(maxEAR > EAR_thresh) %>%
    group_by(site,date) %>%
    filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
    group_by(site,date,ID) %>%
    summarize(EARsum = sum(maxEAR))
  
  yaxis_plot_nums <- plot_dimensions[[i]][1]*j-1
  yaxt <- ifelse(j %in%  (0:4*plot_dimensions[[i]][2] +1),"s","n")
  boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
          log="y",main=paste(chnms),
          las=2,
          ylim=c(1e-5,10),
          cex.axis=axis_text_cex,
          cex.main=title_text_cex,
          yaxt=yaxt)
  mtext(paste(nMixSites,"Sites"),side=3,line=-1,cex=0.7)
  
  
}
mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
mtext(bquote(.(y_label[["y_label"]])),
      side = 2,line=3,outer=TRUE, cex = 0.75)
mtext(side = 1, cex = 0.5,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7"~.(LETTERS[i-1])~":") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label[["y_label"]])  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-chemical mixtures present",
                         "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013 (" * italic(.(y_label$caption)) * ").")),outer=TRUE)

dev.off()
