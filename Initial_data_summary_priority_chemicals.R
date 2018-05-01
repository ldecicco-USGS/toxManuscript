library(dplyr)

source("data_setup.R")

chem_info <- tox_list[[2]]
names(chem_info)[2] <- "chnm"

# Determine which compounds are:
# 1. present at a min of 10 sites
# 2. EAR > thresh at a min of 5 sites

# Set thresholds
thresh <- 0.001 #EAR threshold
num_sites_occur_thresh <- 10
num_sites_exceed_thresh <- 5



#Filter to chemicals that occur at > num_sites_occur_thresh
greater_than_thresh <- filter(chemicalSummary, EAR > thresh) %>%
  group_by(CAS) %>%
  summarize(numSites = n_distinct(shortName),
            EARmax = max(EAR)) %>%
  filter(numSites >= num_sites_exceed_thresh)

#Filter to chemicals that exceed the EAR threshold at > num_sites_exceed_thresh
Num_sites_occur <- group_by(chemicalSummary, CAS) %>%
  filter(EAR > 0) %>%
  summarize(numSitesOccur = n_distinct(site),
            EARmax = max(EAR)) %>%
  filter(numSitesOccur >= num_sites_occur_thresh)

#Combine info for occurrence and exceedance thresholds to develop one data frame
priority_chems <- as.character(Num_sites_occur$CAS[which(Num_sites_occur$CAS %in% greater_than_thresh$CAS)])
priority_chems_df <- Num_sites_occur[which(Num_sites_occur$CAS %in% greater_than_thresh$CAS),]

chem_info_priority_chems <- filter(chem_info,CAS %in% priority_chems) %>%
  arrange(Class,chnm)
chem_info_priority_chems <- left_join(priority_chems_df,chem_info_priority_chems,by="CAS")%>%
arrange(Class,chnm)
write.csv(chem_info_priority_chems,file="priority_chems.csv",row.names = FALSE)

