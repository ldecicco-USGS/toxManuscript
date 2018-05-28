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

ear_thresh <- 0.001
siteThres <- 10
# ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)
relevance$Relevant <- MakeTitles(relevance$Relevant)

relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.) 

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= ear_thresh) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThres)

priority_endpoints <- endpoints_sites_hits$endPoint

boxData_max <- chemicalSummary %>%
  left_join(AOP, by="endPoint") %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)])




######################################################################################
#Code for exploring data to be included in manuscript text
relevantEARs <- filter(boxData, grepl("yes|maybe",Relevant,ignore.case = TRUE)) 
range(relevantEARs$maxMaxEAR) # Get max EAR

# Determine which chemicals for each AOP, range of EARs, and how many sites

#Start with chems identified from Fig1: present at 10 or more sites with EARmaxChem > 10-3 at 5 or more sites
priority_chems <- read.csv("priority_chems.csv",stringsAsFactors = TRUE)

AOPs_by_priority_chem <- filter(boxData_max,CAS %in% priority_chems$CAS) %>%
  left_join(relevance,by=("ID")) %>%
  filter(maxEAR > 0) %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE))%>%
  group_by(ID,chnm,CAS) %>%
  summarize(nSites = n_distinct(site),
            maxEAR = max(maxEAR),
            maxEndpoint = endPoint[which.max(maxEAR)]) %>%
  filter(maxEAR > 0.001) %>%
  arrange(ID,nSites)

unique(AOPs_by_priority_chem$chnm)
priority_chems[which(!priority_chems$CAS %in% AOPs_by_priority_chem$CAS),]


## find most common mixtures for each AOP
mixtures_by_sample <- filter(chemicalSummary, CAS %in% priority_chems$CAS) %>%
  filter(EAR > 0.001) %>%
  left_join(AOP,by="endPoint") %>%
  left_join(relevance,by="ID") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(site,date) %>%
  summarize(chemVector = paste(sort(unique(CAS)),collapse = "; ")) %>%
  group_by(site,chemVector) %>%
  summarize(numOccur = n_distinct(date)) %>%
  arrange(site)

as.data.frame(table(mixtures_by_sample$chemVector))

unique(mixtures_by_sample$chemVector)

## find most common mixtures for each AOP
mixtures_by_sample_AOP <- filter(chemicalSummary, CAS %in% priority_chems$CAS) %>%
  filter(EAR > 0.001) %>%
  left_join(AOP,by="endPoint") %>%
  left_join(relevance,by="ID") %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
  group_by(site,date,ID) %>%
  summarize(chemVector = paste(sort(unique(CAS)),collapse = "; ")) %>%
  group_by(site,ID,chemVector) %>%
  summarize(numOccur = n_distinct(date)) %>%
  arrange(site)

as.data.frame(table(mixtures_by_sample_AOP$ID, mixtures_by_sample_AOP$chemVector))

unique(mixtures_by_sample_AOP$chemVector)

range(test$numOccur)
unique(test$CAS)

####################################################################################
### Mixtures analysis more thoroughly

Chem_vectors_by_site <- filter(chemicalSummary, CAS %in% priority_chems$CAS) %>%
  filter(EAR > 0.0001) %>%
  group_by(site,date) %>%
  summarize(chemVector = paste(sort(unique(CAS)),collapse = "; "))

unique(test$chemVector)  


for(i in 1:dim(priority_chems)[1]) {
  chem <- priority_chems[i,"CAS"]
  sites_by_vector <- filter(Chem_vectors_by_site,grepl(chem,chemVector))
  Num_sites_by_vector <- data.frame(numSites =length(unique(sites_by_vector$site)))
  
  Num_sites_by_vector$chemVector <- chem
  Num_sites_by_vector$nChems <- 1
  if(i==1) Num_sites_by_mixture <- Num_sites_by_vector
  else Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  
  for(j in 2:dim(priority_chems)[1]){
    chem2 <- priority_chems[j,"CAS"]
    if(chem2 != chem){
      chems2 <- paste0(chem,"; ",chem2)
      sites_by_vector <- filter(sites_by_vector,grepl(chem2,chemVector))
      Num_sites_by_vector <- data.frame(numSites =length(unique(sites_by_vector$site)))
      Num_sites_by_vector$chemVector <- chems2
      Num_sites_by_vector$nChems <- 2
      Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
    }
    
    for(k in 3:dim(priority_chems)[1]){
      chem3 <- as.character(priority_chems[k,"CAS"])
      if(!grepl(chem3,chems2)){
        chems3 <- paste0(chems2,"; ",chem3)
        sites_by_vector <- filter(sites_by_vector,grepl(chem3,chemVector))
        Num_sites_by_vector <- data.frame(numSites =length(unique(sites_by_vector$site)))
        Num_sites_by_vector$chemVector <- chems3
        Num_sites_by_vector$nChems <- 3
        Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
      }
      for(l in 3:dim(priority_chems)[1]){
        chem4 <- as.character(priority_chems[l,"CAS"])
        if(!grepl(chem4,chems3)){
          chems4 <- paste0(chems3,"; ",chem4)
          sites_by_vector <- filter(sites_by_vector,grepl(chem4,chemVector))
          Num_sites_by_vector <- data.frame(numSites =length(unique(sites_by_vector$site)))
          Num_sites_by_vector$chemVector <- chems4
          Num_sites_by_vector$nChems <- 4
          Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
        }
      }
    }
  }
}
Num_sites_by_mixture <- filter(Num_sites_by_mixture,numSites>0) %>%
  arrange(nChems,desc(numSites))


#########################################################################################


chemRows <- sum(grepl(chem,Chem_vectors_by_site$chemVector))
unique(Chem_vectors_by_site[chemRows,"site"]
       