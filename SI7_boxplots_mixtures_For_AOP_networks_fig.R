library(toxEval)
library(dplyr)
library(tidyr)
library(data.table)
#########################################################################################

source(file = "Table_SI8_mixtures_v4.R")

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

test1 <- chemSummData_max

EAR_thresh <- 0.00001

#plot_dimensions <- list(c(0,0),c(4,4),c(3,4),c(3,4))
margins <- c(4,1.5,1,0)
outer_margins <- c(7,5,2,1)
axis_text_cex <- 0.6
title_text_cex <- 1
yaxt <- "s"
yaxis <- c(0.00001,0.001,0.1,10)
names(yaxis) <- c("0.00001","0.001","0.1","10.0")

y_label <- bquote(EAR[SiteAOP])

# Add chemical names into ToxCast_ACC df
library(toxEval)
rm(ToxCast_ACC)
ToxCast_ACC <- ToxCast_ACC 
tox_chemicals <- tox_chemicals

ToxCast_ACC <- dplyr::left_join(ToxCast_ACC,
                                dplyr::select(tox_chemicals,
                                              CAS = Substance_CASRN,
                                              chnm = Substance_Name),
                                by="CAS")

ToxCast_ACC$chnm[!is.na(ToxCast_ACC$chnm) & ToxCast_ACC$chnm == "TDCPP"] <- "Tris(1,3-dichloro-2-propyl)phosphate"



###################################
filenm <- "plots/SI7_mixtureBoxplots_AOP_networks.pdf"
pdf(filenm)
par(mfrow=c(3,1),mar=margins,oma=outer_margins)

## Mixture 1 ##
mixture1 <- c("126-73-8","58-08-2")
names(mixture1) <- c("Tributyl phosphate", "Caffeine")
nMixSites <- 14

CASnums <- mixture1
chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])


subChemSummary <- chemSummData_max %>%
  filter(maxEAR > EAR_thresh) %>%
  group_by(site,date) %>%
  filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
  group_by(site,date,ID) %>%
  summarize(EARsum = sum(maxEAR))

length(unique(subChemSummary$site))

test <- subChemSummary %>%
  group_by(site,date) %>%
  summarize(EARsum2 = sum(EARsum)) %>%
  filter(EARsum2 > 0.001)

length(unique(test$site))


AOPs <- sort(unique(subChemSummary$ID))
subChemSummary$ID <- factor(subChemSummary$ID,levels = AOPs,ordered = TRUE)

bp <- boxplot(subChemSummary$EARsum ~ subChemSummary$ID, 
        log="y",
        las=2,
        ylim=c(1e-5,10),
        cex.axis=axis_text_cex,
        cex.main=title_text_cex,
        yaxt="n",xaxt = "n",
        cex.lab=5)
title(paste(chnms, collapse = "\n"), line=-2, cex.main = title_text_cex,adj=0.03)
mtext(paste("EAR > 0.001 at",nMixSites,"Sites"),side=3,line=0,cex=1)

axis(side = 1,at = 1:length(AOPs),labels = AOPs,las=2)
axis(side = 2,at = yaxis,labels = names(yaxis),las=2)
## Mixture 2 ##

mixture2 <- c("80-05-7","134-62-3")
names(mixture2) <- c("Bisphenol A", "DEET")

CASnums <- mixture2
nMixSites <- 19
chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])

subChemSummary <- chemSummData_max %>%
  filter(maxEAR > EAR_thresh) %>%
  group_by(site,date) %>%
  filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
  group_by(site,date,ID) %>%
  summarize(EARsum = sum(maxEAR))

length(unique(subChemSummary$site))

  AOPs <- sort(unique(subChemSummary$ID))
subChemSummary$ID <- factor(subChemSummary$ID,levels = AOPs,ordered = TRUE)

bp <- boxplot(subChemSummary$EARsum ~ subChemSummary$ID, 
              log="y",
              las=2,
              ylim=c(1e-5,10),
              cex.axis=axis_text_cex,
              cex.main=title_text_cex,
              yaxt="n",xaxt = "n",
              cex.lab=5)
title(paste(chnms, collapse = "\n"), line=-2, cex.main = title_text_cex,adj=0.03)
mtext(paste("EAR > 0.001 at",nMixSites,"Sites"),side=3,line=0,cex=1)

axis(side = 1,at = 1:length(AOPs),labels = AOPs,las=2)
axis(side = 2,at = yaxis,labels = names(yaxis),las=2)


## Mixture 3 ##
mixture3 <- c("84852-15-3","1912-24-9", "51218-45-2")
names(mixture3) <- c("4-Nonylphenol, branched", "Atrazine","Metolachlor")

CASnums <- mixture3
nMixSites <- 5
chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% CASnums),"chnm"])

subChemSummary <- chemSummData_max %>%
  filter(maxEAR > EAR_thresh) %>%
  group_by(site,date) %>%
  filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
  group_by(site,date,ID) %>%
  summarize(EARsum = sum(maxEAR))

length(unique(subChemSummary$site))

  AOPs <- sort(unique(subChemSummary$ID))
subChemSummary$ID <- factor(subChemSummary$ID,levels = AOPs,ordered = TRUE)

bp <- boxplot(subChemSummary$EARsum ~ subChemSummary$ID, 
              log="y",
              las=2,
              ylim=c(1e-5,10),
              cex.axis=axis_text_cex,
              cex.main=title_text_cex,
              yaxt="n",xaxt = "n",
              cex.lab=5)
title(paste(chnms, collapse = "\n"), line=-3, cex.main = title_text_cex,adj=0.03)
mtext(paste("EAR > 0.001 at",nMixSites,"Sites"),side=3,line=0,cex=1)


axis(side = 1,at = 1:length(AOPs),labels = AOPs,las=2)
axis(side = 2,at = yaxis,labels = names(yaxis),las=2)


## Add labels and caption to fig

  mtext("AOP ID",side=1,outer=TRUE, line = -1, cex = 1.3)
  # mtext(paste(i,"-Compound Mixtures"),outer=TRUE)
  mtext(bquote(.(y_label)),
        side = 2,line=2.5,outer=TRUE, cex = 1.2)
  mtext(side = 1, cex = 0.5,adj = 0,line = 2,
        text = bquote(atop(bold("Figure SI-7.") ~ "Boxplots of exposure activity ratios (" *
                             .(y_label)  *
                             ") by adverse outcome pathway for three selected chemical mixtures present",
                           "in samples that occurred at a minimum of 4 sites during monitoring of Great Lakes tributaries, 2010-2013.")),outer=TRUE)


dev.off()


