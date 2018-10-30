#Determine increase in EARs due to mixtures as opposed to individual chemicals

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
  group_by(endPoint,site,chnm,CAS,date) %>%
  summarize(EARChem = max(EAR))

  

endpoints_sites_hits2 <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARSiteChem = sum(EAR))

ep_joined <- left_join(endpoints_sites_hits,endpoints_sites_hits2)
ep_joined$EARproportion <- ep_joined$EARChem/ep_joined$EARSiteChem
range(ep_joined$EARproportion)

#unique(ep_joined$endPoint)


#############################################################
# How many chemicals have EARchem>0.001 in the entire data set
ep_chemGT001 <- filter(ep_joined, EARChem > 0.001)
ep_chemGT001$chnm <- as.character(ep_chemGT001$chnm)

unique(ep_chemGT001$chnm)
table(ep_chemGT001$chnm) #23 chemicals

#############################################################
#How many chemicals are present for instances when EARSiteMix > 0.001
ep_endpointGT001 <- filter(ep_joined, EARSiteChem > 0.001) %>%
  group_by(chnm) %>%
  summarize(num_sites = length(unique(site)))
ep_endpointGT001$chnm <- as.character(ep_endpointGT001$chnm)

dim(ep_endpointGT001)[1] #42 chemicals contribute


#Now filter to those chemicals that contribute > X%
#ep_endpointGT001$site <- as.factor(ep_endpointGT001$site)
ep_proportionGT.01 <-  filter(ep_joined, EARSiteChem > 0.001) %>%
  filter(EARproportion > 0.05) %>%
  group_by(chnm) %>%
  summarize(num_sites = length(unique(site)))

dim(ep_proportionGT.01)[1] #38 chemicals contribute
table(ep_proportionGT.01$chnm)


#Now compare the EARsiteMix (summation by endpoint) to the max EARChem
max_chem_contributions_per_sample <-
ep_joined %>% 
  group_by(endPoint,site,date) %>%
  summarize(maxProportion = max(EARproportion)) # EARSiteMix is up to 76% 
                                                # greater than EARChem for the max chem contributor


