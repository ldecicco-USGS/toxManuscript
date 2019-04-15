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
source(file = "combo_graph_function.R")

# Choosing priority chemicals
# ---------------------------
# 1. Join EAR data with AOP relevance information
# 2. Remove AOPs that are not relevant
# 3. Filter by max EAR per chemical per endpoint
# 4. sum EARs by sample and AOP (EARsumAOP)
# 5. Retain only samples with EARsumAOP > 10^-3
# 6. Remove chemicals in each sample that do not contribute at least 1% of EARsumAOP
# 7. Retain only chemicals that show up at a minimum of 5 sites


EAR_thresh <- 0.001
# ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk_Dec_2018.csv", stringsAsFactors = FALSE)
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

#Look at percent range for benzophenone:
range(filter(boxData_pct, chnm == "Benzophenone")$EAR_percent) #result: max =1.4%

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

test <- boxData_pct%>% filter(chnm %in% AOP_priority_chems)
test$chnm <- as.character(test$chnm)

options(scipen=10)
par(mar =c(5,15,2,1))
boxplot(EAR_percent ~ chnm,data=test,log="x",horizontal=TRUE,las=2)
boxplot(EARsum ~ chnm,data=test,log="x",horizontal=TRUE,las=2)

#Check 3,4-dichlorophenyl isocyanate
test <- chemicalSummary %>% filter(CAS == "102-36-3") %>%
  filter(EAR>0.00)
range(test$EAR)
test$chnm <- as.character(test$chnm)
boxplot(EAR ~ chnm,data=test,log="x",horizontal=TRUE,las=2)


####################################################################################
### Thorough mixtures analysis up to 4 chemicals

# Mixtures analysis
# ------------------
# 1. Join EAR data with AOP relevance information
# 2. Remove AOPs that are not relevant
# 3. Filter by max EAR per chemical per endpoint
# 4. sum EARs by sample and AOP (EARsumAOP)
# 5. Retain only samples with EARsumAOP > 10^-3
# 6. Take note of which samples remain
# 7. Go back to original data set and subset only chemicals that result from step 6
# 8. Subset to "priority chemicals" defined above.
# 9. Remove individual instances of EAR < 0.00001 (this is < 1% of potential influence in individual EARsumAOP values) and it makes the resulting data set more manageable for mixture analysis
# 10. Determine how many sites that 2-, 3-, and 4-chemical combinations occur at.
# 11. Examine EARsumAOPs for resulting data

# Using the same code as in fig. 1:
EAR_thresh <- 0.001

# gd <- graph_chem_data(chemicalSummary)
# priority_chems2 <- get_priority_chms(gd, EAR_thresh)
# priority_cas <- tox_chemicals %>%
#   filter(Substance_Name %in% priority_chems2) %>%
#   pull(Substance_CASRN)
#   
# # priority_chems <- read.csv("priority_chems.csv",stringsAsFactors = FALSE) #chems resulting from fig 1 analysis
# AOP_priority_CAS[!AOP_priority_CAS %in% priority_cas]
# # priority_chems[!priority_cas %in% AOP_priority_CAS,"chnm"]

chemSummaryAOP <- boxData_max %>%
  group_by(ID, site, date) %>%
  summarize(EARsum = sum(maxEAR, na.rm = TRUE))%>%
  filter(EARsum > EAR_thresh) %>%
  mutate(sample = paste(site,date))

EAR_thresh_individual_chem <- 0.000001

Chem_vectors_by_site <- chemicalSummary %>%
  mutate(sample = paste(site,date)) %>%
  filter(sample %in% chemSummaryAOP$sample)%>%
  filter(CAS %in% AOP_priority_CAS) %>% 
  filter(EAR > EAR_thresh_individual_chem) %>%
  group_by(site,date) %>%
  summarize(chemVector = paste(sort(unique(CAS)),collapse = "|"))

allSTAIDs1 <- character()
for(i in 1:length(AOP_priority_CAS)) {
  chem <- AOP_priority_CAS[i]
  sites_by_vector <- filter(Chem_vectors_by_site,grepl(chem,chemVector))
  STAIDs <- unique(sites_by_vector$site)
  Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
  Num_sites_by_vector$chemVector <- chem
  Num_sites_by_vector$nChems <- 1
  Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
  allSTAIDs1 <- unique(c(allSTAIDs1,STAIDs))
  if(i==1) {Num_sites_by_mixture <- Num_sites_by_vector
  } else Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  
#  if(i < length(AOP_priority_CAS)){
  allSTAIDs2 <- character()
  for(j in (i):length(AOP_priority_CAS)){
    chem2 <- AOP_priority_CAS[j]
    if(chem2 != chem){
      chems2 <- paste0(sort(AOP_priority_CAS[c(i,j)]),collapse="|")
#      chems2_list <- c(chems2_list,AOP_priority_CAS[c(i,j)]))
      sites_by_vector <- filter(sites_by_vector,grepl(chem2,chemVector)) %>%
        filter(site %in% allSTAIDs1)
      STAIDs <- unique(sites_by_vector$site)
      Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
      Num_sites_by_vector$chemVector <- chems2
      Num_sites_by_vector$nChems <- 2
      Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
      allSTAIDs2 <- unique(c(allSTAIDs2,STAIDs))
      Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
    }
  }
}

chem2_df <- filter(Num_sites_by_mixture,nChems==2 & numSites > 0)


#Filter to 2-chem mixtures with more than zero sites
#loop through 2-chem mixtures to look for 3 chem mixtures

# Determine unique 3-chem combos
chems_char_vector <- character()
for(m in 1:dim(chem2_df)[1]){
  current_chems <- unlist(strsplit(chem2_df$chemVector[m],split = "|",fixed=TRUE))
  other_chems <- AOP_priority_CAS[-which(AOP_priority_CAS %in% current_chems)]
  for(l in 1:length(other_chems)) {
    chems <- c(current_chems,other_chems[l])
    chems_char_vector <- c(chems_char_vector,paste(sort(c(chems)),collapse="|"))
  }
}
length(chems_char_vector)

chems_char_vector <- unique(chems_char_vector)
length(chems_char_vector)

for(m in 1:length(chems_char_vector)){
  chems_char <- chems_char_vector[m]
  chems <- unlist(strsplit(chems_char,split = "|",fixed=TRUE))
  mixture_df <- Chem_vectors_by_site
  for (i in  1:3){
    mixture_df <- mixture_df[grep(chems[i], mixture_df$chemVector),]
  }
  # mixture_rows <- intersect(intersect(
  #   grep(chems[1],Chem_vectors_by_site$chemVector),
  #   grep(chems[2],Chem_vectors_by_site$chemVector)),
  #   grep(chems[3],Chem_vectors_by_site$chemVector))
  #mixture_df <- Chem_vectors_by_site[mixture_rows,]
  unique(mixture_df$site)
  STAIDs <- unique(mixture_df$site)
  Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
  Num_sites_by_vector$chemVector <- chems_char
  Num_sites_by_vector$nChems <- 3
  Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
  
  Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  
}






    # allSTAIDs3 <- character()
    # for(k in j:length(AOP_priority_CAS)){
    #   chem3 <- AOP_priority_CAS[k]
    #   if(all(!duplicated(AOP_priority_CAS[c(i,j,k)]))){
    #     chems3 <- paste0(sort(AOP_priority_CAS[c(i,j,k)]),collapse="|")
    #     sites_by_vector <- filter(sites_by_vector,grepl(chem3,chemVector)) %>%
    #       filter(site %in% allSTAIDs2)
    #     STAIDs <- unique(sites_by_vector$site)
    #     Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
    #     Num_sites_by_vector$chemVector <- chems3
    #     Num_sites_by_vector$nChems <- 3
    #     Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
    #     allSTAIDs3 <- unique(c(allSTAIDs3,STAIDs))
    #     Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
    #   }

  #     # allSTAIDs4 <- character()
  #     # for(l in k:length(AOP_priority_CAS)){
  #     #   chem4 <- AOP_priority_CAS[l]
  #     #   if(all(!duplicated(AOP_priority_CAS[c(i,j,k,l)]))){
  #     #     chems4 <- paste0(sort(AOP_priority_CAS[c(i,j,k,l)]),collapse="|")
  #     #     sites_by_vector <- filter(sites_by_vector,grepl(chem4,chemVector)) %>%
  #     #       filter(site %in% allSTAIDs3)
  #     #     STAIDs <- unique(sites_by_vector$site)
  #     #     Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
  #     #     Num_sites_by_vector$chemVector <- chems4
  #     #     Num_sites_by_vector$nChems <- 4
  #     #     Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
  #     #     allSTAIDs4 <- unique(c(allSTAIDs4,STAIDs))
  #     #     Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  #     #   }
  #     #   allSTAIDs5 <- character()
  #     #   for(m in l:length(AOP_priority_CAS)){
  #     #     chem5 <- AOP_priority_CAS[m]
  #     #     if(all(!duplicated(AOP_priority_CAS[c(i,j,k,l,m)]))){
  #     #       chems5 <- paste0(sort(AOP_priority_CAS[c(i,j,k,l,m)]),collapse="|")
  #     #       sites_by_vector <- filter(sites_by_vector,grepl(chem5,chemVector)) %>%
  #     #         filter(site %in% allSTAIDs4)
  #     #       STAIDs <- unique(sites_by_vector$site)
  #     #       Num_sites_by_vector <- data.frame(numSites =length(STAIDs))
  #     #       Num_sites_by_vector$chemVector <- chems5
  #     #       Num_sites_by_vector$nChems <- 5
  #     #       Num_sites_by_vector$STAIDs <- paste(STAIDs,collapse = "|")
  #     #       allSTAIDs5 <- unique(c(allSTAIDs5,STAIDs))
  #     #       Num_sites_by_mixture <- rbind(Num_sites_by_mixture,Num_sites_by_vector)
  #     #     }
  #     #   }
  #     # }
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
                  paste(siteList[strsplit(Num_sites_by_mixture$STAIDs[i],"\\|")[[1]]],collapse="|"))
}
Num_sites_by_mixture$siteVector <- siteColumn

####Add chemical names from CAS vector
chemList <- as.character(tox_list[["chem_info"]]$"Chemical Name")
names(chemList) <- tox_list[["chem_info"]]$CAS

chemColumn <- character()
for(i in 1:dim(Num_sites_by_mixture)[1]){
  chemColumn <- c(chemColumn,
                  paste(chemList[strsplit(Num_sites_by_mixture$chemVector[i],"\\|")[[1]]],collapse="|"))
}
Num_sites_by_mixture$chnmVector <- chemColumn


write.csv(Num_sites_by_mixture,file="SI_table7 Num_sites_by_mixture_temp.csv",row.names = FALSE)

TableSI8 <- Num_sites_by_mixture[,c("nChems","numSites","chemVector","chnmVector","siteVector","STAIDs")]
names(TableSI8) <- c("Number of Chemicals","Number of sites", "CAS#","Chemical Names","Site Short Names","USGS Station IDs")

write.csv(TableSI8,file="tables/SI_table8_Num_sites_by_mixture.csv",row.names = FALSE)
