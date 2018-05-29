library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(dataRetrieval)
library(cowplot)
library(grid)

####################################
source(file = "data_setup.R")
source(file = "MakeTitles.R")

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


boxData_max <- chemicalSummary %>%
  filter(EAR > 0) %>%
  left_join(AOP_relevance, by="endPoint") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)]) %>%
  mutate(sample = paste(site,date))

boxData_pct <- boxData_max %>%
  group_by(ID,sample) %>%
  summarise(EARsum = sum(maxEAR)) %>%
  right_join(boxData_max,by=c("ID","sample")) %>%
  filter(EARsum > EAR_thresh) %>%
  mutate(EAR_percent = maxEAR/EARsum)

pct_check <- boxData_pct %>%
  group_by(ID,sample) %>%
  summarise(pct_sum = sum(EAR_percent))
range(pct_check$pct_sum) # All good. pct_sum = 1.0

siteThresh <- 5
pct_thresh <- 0.01
filtered_chems <- boxData_pct %>%
  filter(EAR_percent > pct_thresh) %>%
  group_by(chnm,CAS) %>%
  summarize(nSites = n_distinct(site)) %>%
  filter(nSites > siteThresh)

unique(filtered_chems$chnm)

AOP_priority_chems <- as.character(unique(filtered_chems$chnm))
AOP_priority_CAS <-unique(filtered_chems$CAS)

boxplot(EAR_percent ~ chnm,data=boxData_pct,log="x",horizontal=TRUE,las=2)


####################################################################################
### Thorough mixtures analysis up to 5 chemicals
priority_chems <- read.csv("priority_chems.csv",stringsAsFactors = FALSE)
AOP_priority_CAS[!AOP_priority_CAS %in% priority_chems$CAS]
priority_chems[!priority_chems$CAS %in% AOP_priority_CAS,"chnm"]

EAR_thresh <- 0.001
chemSummaryAOP <- chemicalSummary %>%
  filter(EAR > 0) %>%
  left_join(AOP_relevance, by="endPoint") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(ID, site, date) %>%
  summarize(EARsum = sum(EAR, na.rm = TRUE))%>%
  filter(EARsum > EAR_thresh) %>%
  mutate(sample = paste(site,date))

EAR_thresh_individual_chem <- 0.0001
Chem_vectors_by_site <- chemicalSummary %>%
  mutate(sample = paste(site,date)) %>%
  filter(sample %in% chemSummaryAOP$sample)%>%
  filter(CAS %in% AOP_priority_CAS) %>% 
  filter(EAR > 0.000001) %>%
  group_by(site,date) %>%
  summarize(chemVector = paste(sort(unique(CAS)),collapse = "; "))

allSTAIDs1 <- character()
for(i in 1:length(AOP_priority_CAS)) {
  chem <- AOP_priority_CAS[i]
  sites_by_vector <- filter(Chem_vectors_by_site,grepl(chem,chemVector))
  STAIDs <- unique(sites_by_vector$site)
  Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
  Num_sites_by_vector$chemVector <- chem
  Num_sites_by_vector$nChems <- 1
  Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "; ")
  allSTAIDs1 <- unique(c(allSTAIDs1,STAIDs))
  if(i==1) Num_sites_by_mixture <- Num_sites_by_vector
  else Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  
  allSTAIDs2 <- character()
  for(j in i:length(AOP_priority_CAS)){
    chem2 <- AOP_priority_CAS[j]
    if(chem2 != chem){
      chems2 <- paste0(sort(AOP_priority_CAS[c(i,j)]),collapse="; ")
      sites_by_vector <- filter(sites_by_vector,grepl(chem2,chemVector)) %>%
        filter(site %in% allSTAIDs1)
      STAIDs <- unique(sites_by_vector$site)
      Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
      Num_sites_by_vector$chemVector <- chems2
      Num_sites_by_vector$nChems <- 2
      Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "; ")
      allSTAIDs2 <- unique(c(allSTAIDs2,STAIDs))
      Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
    }
    
    allSTAIDs3 <- character()
    for(k in i:length(AOP_priority_CAS)){
      chem3 <- AOP_priority_CAS[k]
      if(all(!duplicated(AOP_priority_CAS[c(i,j,k)]))){
      chems3 <- paste0(sort(AOP_priority_CAS[c(i,j,k)]),collapse="; ")
        sites_by_vector <- filter(sites_by_vector,grepl(chem3,chemVector)) %>%
          filter(site %in% allSTAIDs2)
        STAIDs <- unique(sites_by_vector$site)
        Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
        Num_sites_by_vector$chemVector <- chems3
        Num_sites_by_vector$nChems <- 3
        Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "; ")
        allSTAIDs3 <- unique(c(allSTAIDs3,STAIDs))
        Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
      }
      allSTAIDs4 <- character()
      for(l in i:length(AOP_priority_CAS)){
        chem4 <- AOP_priority_CAS[l]
        if(all(!duplicated(AOP_priority_CAS[c(i,j,k,l)]))){
          chems4 <- paste0(sort(AOP_priority_CAS[c(i,j,k,l)]),collapse="; ")
          sites_by_vector <- filter(sites_by_vector,grepl(chem4,chemVector)) %>%
            filter(site %in% allSTAIDs3)
          STAIDs <- unique(sites_by_vector$site)
          Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
          Num_sites_by_vector$chemVector <- chems4
          Num_sites_by_vector$nChems <- 4
          Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "; ")
          allSTAIDs4 <- unique(c(allSTAIDs4,STAIDs))
          Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
        }
      }
    }
  }
}
        
        
Num_sites_by_mixture <- filter(Num_sites_by_mixture,numSites>0) %>%
  group_by(chemVector,STAIDs) %>%
  summarize(numSites =max(numSites),
            nChems = max(nChems)) %>%
  arrange(nChems,desc(numSites))


####Add site names from siteID vector
siteList <- as.character(tox_list[["chem_site"]]$"Short Name")
names(siteList) <- tox_list[["chem_site"]]$SiteID

siteColumn <- character()
for(i in 1:dim(Num_sites_by_mixture)[1]){
  siteColumn <- c(siteColumn,
                  paste(siteList[strsplit(Num_sites_by_mixture$STAIDs[i],"; ")[[1]]],collapse="; "))
}
Num_sites_by_mixture$siteVector <- siteColumn

####Add chemical names from CAS vector
chemList <- as.character(tox_list[["chem_info"]]$"Chemical Name")
names(chemList) <- tox_list[["chem_info"]]$CAS

chemColumn <- character()
for(i in 1:dim(Num_sites_by_mixture)[1]){
  chemColumn <- c(chemColumn,
                  paste(chemList[strsplit(Num_sites_by_mixture$chemVector[i],"; ")[[1]]],collapse="; "))
}
Num_sites_by_mixture$chnmVector <- chemColumn


#write.csv(Num_sites_by_mixture,file="Num_sites_by_mixture.csv",row.names = FALSE)

#########################################################################################

Num_sites_by_mixture <- read.csv(file="Num_sites_by_mixture.csv",stringsAsFactors = FALSE)
test <- filter(chemicalSummary,shortName=="BlackOH") %>%
  group_by(site,date) %>%
  summarize(EARsum = sum(EAR))


#find samples with the chemicals in each defined mixture
#Add AOP info
#run boxplots for each row in Num_sites_by_mixture. one page for 2 chems, one page for 3...

Chem_vectors_by_site <- filter(chemicalSummary, CAS %in% AOP_priority_CAS) %>% #priority_chems$CAS) %>%
  filter(EAR > 0.000001) %>%
  group_by(site,date) %>%
  summarize(chems = paste0(unique(chnm),collapse=";"))

for(i in 2:5){
  sub_Num_sites <- Num_sites_by_mixture %>%
    filter(numSites == i)
  for(j in 1:dim(sub_Num_sites)[1]){
    CASnums <- strsplit(sub_Num_sites$CAS,"; ")
    subChemSummary <- chemicalSummary %>%
      group_by(site,date) %>%
      filter(all(CASnums %in% CAS))
    
    