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

AOP_crosswalk <- read.csv("AOP_crosswalk_Dec_2018.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)
relevance$Relevant <- MakeTitles(relevance$Relevant)

relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.) 

# Determine max EAR per chemical per endpoint per sample
endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,chnm,CAS,date) %>%
  summarize(EARChem = max(EAR))

  
# Determine summation of EARs per site per endpoint per sample
endpoints_sites_hits2 <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARSiteChem = sum(EAR))

# Compute percent contribution to EARs by endpoint for each chemical
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

#How many individual samples have EAR > 0.001
unique_samples <- paste(ep_chemGT001$site,ep_chemGT001$date)
length(unique(unique_samples)) #551 samples with EAR > 0.001

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

dim(ep_proportionGT.01)[1] #30 chemicals contribute at least 5%
table(ep_proportionGT.01$chnm)


# Compare the EARsiteMix (summation by endpoint) to the max EARChem
max_chem_contributions_per_sample <-
ep_joined %>% 
  group_by(endPoint,site,date) %>%
  summarize(maxProportion = max(EARproportion)) # EARSiteMix is up to 76% 
                                                # greater than EARChem for the max chem contributor
summary(max_chem_contributions_per_sample$maxProportion)

# Compare the EARsiteAOP (summation by AOP) to the max EARChem
max_chem_contributions_per_sample <-
  ep_joined %>% 
  group_by(endPoint,site,date) %>%
  summarize(maxProportion = max(EARproportion)) # EARSiteMix is up to 76% 
# greater than EARChem for the max chem contributor


#Determine how many endpoints have EAR > 10^-3 at more than 10 sites from individual chemicals and
#then determine how many endpoints have EAR > 10^-3 at more than 10 sites with chemical mixtures

# Individual chemicals
# 1. Filter to EAR > 10^-3 for individual chemicals
# 2. Count unique endpoints
# 3. Add AOPs by endpoint
# 4. Count unique AOPs
# Results: 44 unique endPoints
# Results: 50 unique AOPs

endpoint_site_count_chemical <- endpoints_sites_hits %>%
  filter(EARChem>0.001) %>%
  group_by(endPoint) %>%
  summarize(Num_sites = length(unique(site))) %>%
  filter(Num_sites >= 10)

AOP_site_count_chemical <- left_join(endpoint_site_count_chemical,AOP,by = c("endPoint","endPoint"))

length(unique(AOP_site_count_chemical$endPoint))
length(unique(AOP_site_count_chemical$ID))

# #without site filter
# chemSumEARGTThresh <- chemicalSummary %>%
#   filter(EAR >= 0.001)
# 
# unique(chemSumEARGTThresh$endPoint) #144 assays
# unique(chemSumEARGTThresh$chnm) #23 chemicals


# Chemical mixtures
# 1. Filter to EAR > 10^-3 for EARSiteMix
# 2. Count unique endpoints
# 3. Add AOPs by endpoint
# 4. Count unique AOPs
# Results: 46 unique endPoints
# Results: 50 unique AOPs

endpoint_site_count_mixture <- endpoints_sites_hits2 %>%
  filter(EARSiteChem > 0.001) %>%
  group_by(endPoint) %>%
  summarize(Num_sites =  length(unique(site))) %>%
  filter (Num_sites >= 10)

AOP_site_count_mixture <- left_join(endpoint_site_count_mixture,AOP,by = c("endPoint","endPoint"))

length(unique(AOP_site_count_mixture$endPoint))
length(unique(AOP_site_count_mixture$ID))

# #Without site filter
# chemSumEARAssayGTThresh <- chemicalSummary %>%
#   group_by(site,date,endPoint) %>%
#   summarize(
#     EARsum = sum(EAR),
#     chems = paste(unique(chnm),collapse="SPLITME")) %>%
#   filter(EARsum >= 0.001)
# chemSumEARAssayGTThresh <- as.data.frame(chemSumEARAssayGTThresh)
# 
# mixtureChems <- character()
# for(i in 1:dim(chemSumEARAssayGTThresh)[1]){
#   IndAssayChems <- strsplit(chemSumEARAssayGTThresh[i,"chems"],split = "SPLITME")
#   mixtureChems <- unique(c(mixtureChems,IndAssayChems[[1]]))
# }
#                          
# mixtureChems #47 chemicals
# unique(chemSumEARAssayGTThresh$endPoint) #146 Assays

## AOPs EAR comparisons to individual chemicals##
# Determine max EAR per chemical per AOP per sample

chemicalSummary_AOPs <- left_join(chemicalSummary,AOP,by=c("endPoint","endPoint"))

AOPs_sites_hits_max_chm_by_endpoint <- 
  filter(chemicalSummary_AOPs,EAR > 0) %>%
  group_by(ID,site,chnm,CAS,date) %>%
  summarize(EARChem = max(EAR)) 

AOPs_sites_hits <- AOPs_sites_hits_max_chm_by_endpoint %>%
  group_by(ID,site,date) %>%
  summarize(EARSiteAOP = sum(EARChem))

# Join EARChem info with EARSiteAOP info and determine proportion contribution 
# from individual chemicals
AOPs_site_hits2 <- left_join(AOPs_sites_hits_max_chm_by_endpoint,AOPs_sites_hits) %>%
  filter(EARSiteAOP > 0.001)

test <-              #51 unique AOPS is consistent with previous counts
AOPs_site_hits2 %>%
  group_by(ID) %>%
  summarize(NumSites = length(unique(site))) %>%
  filter(NumSites >=10)

# test to make sure we are accounting for all EARs: all of the EARSumSum values shoudl be 1
AOPs_chem_comparison_check <-
AOPs_site_hits2 %>% 
  group_by(ID,site,date) %>%
  summarize(EARsum = sum(EARChem),
            EARAOPmax = max(EARSiteAOP))
AOPs_chem_comparison_check$EARSumSum <- AOPs_chem_comparison_check$EARsum/AOPs_chem_comparison_check$EARAOPmax
range(AOPs_chem_comparison_check$EARSumSum)

#Determine max chem for each AOP, and then determine how much EAR increase there is
#for the AOP chemical mixture
## Result: EARs for AOPs increase as much as 73% for EARSiteAOP over the max EARChem per AOP
AOPs_chem_comparison <- 
AOPs_site_hits2 %>% 
  group_by(ID,site,date) %>%
  summarize(EARChemMax = max(EARChem),
            EARSiteAOP = max(EARSiteAOP)) %>%
  filter(EARChemMax > 0.001)


AOPs_chem_comparison$EARProportion <- AOPs_chem_comparison$EARChemMax/AOPs_chem_comparison$EARSiteAOP

EARSiteAOPincrease <- (AOPs_chem_comparison$EARSiteAOP - AOPs_chem_comparison$EARChemMax)/AOPs_chem_comparison$EARSiteAOP

summary(EARSiteAOPincrease)
summary(AOPs_chem_comparison$EARProportion)

summary(ep_joined$EARproportion) # 

summary(AOPs_chem_comparison$EARChemMax)
summary(AOPs_chem_comparison$EARSiteAOP)
summary(AOPs_chem_comparison$EARProportion)  

range(AOPs_chem_comparison$EARProportion)  
range(AOPs_chem_comparison$EARSiteAOP)
plot(AOPs_chem_comparison$EARProportion)


#Now filter to those chemicals that contribute > X% for EARSiteAOP
#ep_endpointGT001$site <- as.factor(ep_endpointGT001$site)
AOPs_site_hits2$EARproportion <- AOPs_site_hits2$EARChem/AOPs_site_hits2$EARSiteAOP
ep_proportionGT.01 <-  filter(AOPs_site_hits2,EARSiteAOP > 0.001) %>%
  filter(EARproportion > 0.05) %>%
  group_by(chnm) %>%
  summarize(num_sites = length(unique(site)))

dim(ep_proportionGT.01)[1] #38 chemicals contribute
table(ep_proportionGT.01$chnm)



# What is the contribution of individual chemicals to EARSiteMix (summed by endpoint)

# Determine max EAR per chemical per endpoint per sample
endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,chnm,CAS,date) %>%
  summarize(EARChem = max(EAR))

# Determine summation of EARs per site per endpoint per sample
endpoints_sites_hits2 <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARSiteChem = sum(EAR))

# Compute percent contribution to EARs by endpoint for each chemical
ep_joined <- left_join(endpoints_sites_hits,endpoints_sites_hits2)
ep_joined$EARproportion <- ep_joined$EARChem/ep_joined$EARSiteChem
range(ep_joined$EARproportion)
summary(ep_joined$EARproportion)

# What is the contribution of individual chemicals to EARAOPMix (summed by AOP)

chemicalSummary_AOP <- left_join(chemicalSummary,AOP,by = c("endPoint","endPoint"))

# Determine max EAR per chemical per endpoint per sample
AOP_sites_hits <- filter(chemicalSummary_AOP,EAR > 0) %>%
  group_by(ID,site,chnm,CAS,date) %>%
  summarize(EARChem = max(EAR))

# Determine summation of EARs per site per endpoint per sample
AOP_sites_hits2 <- filter(chemicalSummary_AOP,EAR > 0) %>%
  group_by(ID,site,date) %>%
  summarize(EARAOPChem = sum(EAR))

# Compute percent contribution to EARs by endpoint for each chemical
AOP_joined <- left_join(AOP_sites_hits,AOP_sites_hits2)
AOP_joined$EARproportion <- AOP_joined$EARChem/AOP_joined$EARAOPChem
range(AOP_joined$EARproportion)
summary(AOP_joined$EARproportion)

unique(AOP_sites_hits[which(AOP_sites_hits$EARChem > 0.001),"chnm"])

EARproportionGTpctThresh <- filter(AOP_joined,EARproportion>0.05)
unique(EARproportionGTpctThresh[which(EARproportionGTpctThresh$EARAOPChem > 0.001),"chnm"])

