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


Num_sites_by_mixture <- read.csv(file="Num_sites_by_mixture.csv",stringsAsFactors = FALSE)

chemSummData_max <- chemicalSummary %>%
  filter(EAR > 0) %>%
  left_join(AOP_relevance, by="endPoint") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE)) 

filenm <- "mixtureBoxplotsChosen.pdf"
pdf(filenm)
EAR_thresh <- 0.00001


plot_dimensions <- c(2,3)
margins <- c(2,0.5,4,0)
outer_margins <- c(2,4,2,1)
axis_text_cex <- 0.6
title_text_cex <- 0.7

par(mfrow=plot_dimensions,mar=margins,oma=outer_margins)
plot_count <- 0

stat_stuff <- list()
for(i in 2:3) {
  sub_Num_sites <- Num_sites_by_mixture %>%
    filter(nChems == i,numSites>=4)
  for(j in 1:dim(sub_Num_sites)[1]){
    plot_condition <- (i==2 & j %in% c(3,5)) | (i ==3 & j %in% c(2,5,6))
    if(plot_condition) {
      plot_count <- plot_count + 1
      CASnums <- strsplit(sub_Num_sites$chemVector[j],"; ")[[1]]
      nMixSites <- sub_Num_sites[j,"numSites"]
      chnms <- unique(as.data.frame(ACC)[which(ACC$casn %in% CASnums),"chnm"])
      
      subChemSummary <- chemSummData_max %>%
        filter(maxEAR > EAR_thresh) %>%
        group_by(site,date) %>%
        filter(grepl(paste(CASnums,collapse="|"), CAS)) %>%
        group_by(site,date,ID) %>%
        summarize(EARsum = sum(maxEAR))
      
      yaxt <- ifelse(plot_count %in% c(1,4),"s","n")
      
      mean_graph <- subChemSummary %>%
        mutate(ID = as.character(ID)) %>%
        group_by(ID) %>%
        summarize(mean = mean(EARsum, na.rm = TRUE))
      
      bp_val <- boxplot(subChemSummary$EARsum ~ as.character(subChemSummary$ID), 
              log="y",main=paste(chnms),
              las=2,
              ylim=c(1e-5,10),
              cex.axis=axis_text_cex,
              cex.main=title_text_cex,
              yaxt=yaxt)
      

      stat_df <- data.frame(median = bp_val$stats[3,],
                            max = bp_val$stats[5,], 
                            ID = bp_val$names,
                            stringsAsFactors = FALSE) %>%
        left_join(mean_graph, by="ID") %>%
        select(ID, median, mean, max)
      
      stat_stuff[[paste(chnms, collapse = ";")]] <- stat_df

      mtext(paste(nMixSites,"Sites"),side=3,line=-1,cex=0.7)
      mtext("AOP",side=1,outer=TRUE)
      mtext("EAR sum by sample",side = 2,line=2.5,outer=TRUE)
      if(plot_count==1){
        mixtureAOPs = data.frame(paste(chnms,collapse = "; "),
                                 paste(unique(as.character(subChemSummary$ID)),collapse="; "))
      } else {
        mixtureAOPs = rbind(mixtureAOPs,data.frame(paste(chnms,collapse = "; "),
                                                   paste(unique(as.character(subChemSummary$ID)),collapse="; ")))
      }
    }
  }
  mtext(paste("Chosen Mixtures for AOP Network Development"),outer=TRUE)
}



dev.off()
shell.exec(filenm)
names(mixtureAOPs) <- c("chnms","AOP_IDs")
write.csv(mixtureAOPs, file="AOP_IDs_for_priority_mixtures.csv",row.names = FALSE)

library(openxlsx)
write.xlsx(stat_stuff, file = "EAR_stats.xlsx", append=TRUE)
