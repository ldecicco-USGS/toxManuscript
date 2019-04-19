library(toxEval)
library(dplyr)
library(tidyr)
library(data.table)
#########################################################################################
source(file = "data_setup.R")
#source(file = "MakeTitles.R")
#source(file = "combo_graph_function.R")

mixtures <- read.csv("SI_table7 Num_sites_by_mixture_temp.csv",stringsAsFactors = FALSE)
#source(file = "Table_SI8_mixtures_v4.R")

#Initialize some things
#plot_dimensions <- list(c(0,0),c(4,4),c(3,4),c(3,4))
margins <- c(4,1.5,1,0)
outer_margins <- c(7,5,2,1)
axis_text_cex <- 0.6
title_text_cex <- 1
yaxt <- "s"
yaxis <- c(0.00001,0.001,0.1,10)
names(yaxis) <- c("0.00001","0.001","0.1","10.0")
textline <- c(-2,-2,-3)

y_label <- bquote(EAR[SiteAOP])

# Add chemical names into ToxCast_ACC df
rm(ToxCast_ACC)
ToxCast_ACC <- ToxCast_ACC 
tox_chemicals <- tox_chemicals

ToxCast_ACC <- dplyr::left_join(ToxCast_ACC,
                                dplyr::select(tox_chemicals,
                                              CAS = Substance_CASRN,
                                              chnm = Substance_Name),
                                by="CAS")

ToxCast_ACC$chnm[!is.na(ToxCast_ACC$chnm) & ToxCast_ACC$chnm == "TDCPP"] <- "Tris(1,3-dichloro-2-propyl)phosphate"



#1. Compute EARAOP values for each sample.
#2. For each of the three mixtures, filter to just data from mixtures
#3. Generate boxplots for each of the three mixtures

# EAR_thresh <- 0.001
mixture_CAS <- list(c("126-73-8","58-08-2"),c("80-05-7","134-62-3"),c("84852-15-3","1912-24-9", "51218-45-2"))
names(mixture_CAS[[1]]) <- c("Tributyl phosphate", "Caffeine")
names(mixture_CAS[[2]]) <- c("Bisphenol A", "DEET")
names(mixture_CAS[[3]]) <- c("4-Nonylphenol, branched", "Atrazine","Metolachlor")

#Sum all EARAOPs by sample and include a chemvector and EARsum
#Use this df for determining how many sites a mixture is present regardless of EAR level

###################################
# Make boxplots                   #
###################################
filenm <- "plots/SI7_mixtureBoxplots_AOP_networks.pdf"
#filenm <- "SI7_mixtureBoxplots_AOP_networks.pdf"
pdf(filenm)
par(mfrow=c(3,1),mar=margins,oma=outer_margins)

for (i in 1:length(mixture_CAS)){
  chems <- mixture_CAS[[i]]
  
  chemPresentAllSamples <- chemicalSummary %>%
    filter(EAR>0,CAS %in% chems)%>%
    left_join(AOP_relevance, by="endPoint") %>%
    filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
    group_by(ID, chnm, CAS, site, date) %>%
    summarize(maxEAR = max(EAR, na.rm = TRUE),
              endPoint_used = endPoint[which.max(EAR)])
  
  chemInSamples <- chemPresentAllSamples %>%
    group_by(site,date) %>%
    summarize(nChems = length(unique(CAS))) %>%
    filter(nChems == length(chems)) %>%
    mutate(sample = paste(site,date))
  samples <- chemInSamples$sample
  
  mixtureEARs <- chemPresentAllSamples %>%
    mutate(sample = paste(site,date)) %>%
    filter(sample %in% samples) %>%
    group_by(ID,site,date,sample) %>%
    summarize(EARsum = sum(maxEAR, na.rm = TRUE))
  
  chnms <- unique(as.data.frame(ToxCast_ACC)[which(ToxCast_ACC$CAS %in% chems),"chnm"])
  
  nMixSites <- length(unique(mixtureEARs$site))
  
  AOPs <- sort(unique(mixtureEARs$ID))
  mixtureEARs$ID <- factor(mixtureEARs$ID,levels = AOPs,ordered = TRUE)
  
  bp <- boxplot(mixtureEARs$EARsum ~ mixtureEARs$ID, 
                log="y",
                las=2,
                ylim=c(1e-5,10),
                cex.axis=axis_text_cex,
                cex.main=title_text_cex,
                yaxt="n",xaxt = "n",
                cex.lab=5)
  title(paste(chnms, collapse = "\n"), line=textline[i], cex.main = title_text_cex,adj=0.03)
  mtext(paste("Mixture present at",nMixSites,"Sites"),side=3,line=0,cex=1)
  
  axis(side = 1,at = 1:length(AOPs),labels = AOPs,las=2)
  axis(side = 2,at = yaxis,labels = names(yaxis),las=2)
  
}


mtext("AOP ID",side=1,outer=TRUE, line = -1.5, cex = 0.65)
# mtext(paste(i,"-Compound Mixtures"),outer=TRUE)
mtext(bquote(paste('EAR'['SiteAOP'])),side=2,line=3,outer=TRUE)
mtext(side = 1, cex = 0.75,adj = 0,line = 2,
      text = bquote(atop(bold("Figure SI-7:") ~ "Boxplots of exposure activity ratios (" *
                           .(y_label)  *
                           ") by adverse outcome pathway for" ~ .(i) ~ "-",
                         )),outer=TRUE)
mtext("chemical mixtures present in samples of Great Lakes tributaries, 2010-2013.",
      outer=TRUE,side=1,line=3,cex=0.75,adj=0)
dev.off()
#shell.exec(filenm)

