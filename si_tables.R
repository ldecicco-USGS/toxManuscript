library(toxEval)
library(dplyr)
library(tidyr)
library(readxl)
library(data.table)

#########################################
## SI 1&2

# si_1 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-1 Watershed Characteristics")
# si_2 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-2 Chemical Classes")

####################################
source(file = "data_setup.R")

cas <- readxl::read_excel("Supplemental_v16.xlsx", sheet = "SI-2 Chemical Classes", skip = 4) %>%
  select(CAS = `...3`) %>%
  filter(!is.na(CAS))
  
cas_key <- toxEval::tox_chemicals %>%
  select(CAS = Substance_CASRN, chnm_tox = Substance_Name) %>%
  full_join(select(tox_list$chem_info, CAS, orig_name = `Chemical Name`), by="CAS") %>%
  filter(CAS %in% tox_list$chem_info$CAS) 

cas_key$chnm_tox[is.na(cas_key$chnm_tox)] <- cas_key$orig_name[is.na(cas_key$chnm_tox)]
cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol monoethoxylate, (sum of all isomers; NP1EO)"] <- paste0("M",cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol monoethoxylate, (sum of all isomers; NP1EO)"])
cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol diethoxylate  (sum of all isomers; NP2EO)"] <- paste0("M",cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol diethoxylate  (sum of all isomers; NP2EO)"])
cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol diethoxylate (OP2EO)"] <- paste0("M",cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol diethoxylate (OP2EO)"])
cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol monoethoxylate (OP1EO)"] <- paste0("M",cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol monoethoxylate (OP1EO)"])
cas_key$chnm <- cas_key$chnm_tox

cas_key_ordered <- cas %>%
  left_join(select(cas_key, CAS, chnm), by="CAS")

dir.create("tables", showWarnings = FALSE)
write.csv(cas_key_ordered, file = "tables/cas_names.csv", row.names = FALSE, na = "")

#########################################
# SI 3:
# 
endPointInfo <- clean_endPoint_info(end_point_info)

si_3_endpoints <- select(chemicalSummary, 
                         `ToxCast Endpoint` = endPoint) %>%
  distinct() %>%
  left_join(select(endPointInfo, 
                   `ToxCast Endpoint`=assay_component_endpoint_name, 
                   `Assay Source`=assay_source_long_name),by = "ToxCast Endpoint") %>%
  distinct() 

# si_3_endpoints <- si_3_endpoints %>%
#   left_join(select(AOP_crosswalk, 
#                    endPoint=`Component Endpoint Name`, AOP = `AOP #`), by=c("ToxCast Endpoint"="endPoint"))
# 
# si_3_endpoints <- si_3_endpoints %>%
#   left_join(select(AOP_info, AOP, Relevant), by="AOP")

dir.create("tables", showWarnings = FALSE)
write.csv(si_3_endpoints, file = "tables/SI3.csv", row.names = FALSE, na = "")
rm(list=ls())
#########################################
# SI 4:

# si_table_counts.R

#########################################
# SI 5:
#
AOP_crosswalk <- fread("AOP_crosswalk_Dec_2018.csv")

write.csv(AOP_crosswalk, file = "tables/SI5.csv", row.names = FALSE, na = "")
rm(list=ls())
#########################################
# SI 6:
#
library(toxEval)
library(dplyr)
library(tidyr)
library(readxl)

source(file = "data_setup.R")
AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")

ear_thresh <- 0.001
siteThres <- 10
# ep_percent_thres <- 0.5

AOP <- read.csv("AOP_crosswalk_Dec_2018.csv", stringsAsFactors = FALSE) %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- data.table::fread("AOP_relevance.csv", data.table = FALSE) %>%
  distinct()

eps_with_ids <- unique(AOP$endPoint)

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= ear_thresh) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThres) %>%
  mutate(hasAOP = endPoint %in% eps_with_ids)

priority_endpoints <- endpoints_sites_hits$endPoint[endpoints_sites_hits$hasAOP]

AOP_full_info <- relevance %>%
  select(-`Endpoint(s)`) %>%
  left_join(AOP, by=c("AOP"="ID")) %>%
  left_join(select(AOP_info, `Abbreviated AOP description` = `...5`, AOP), by="AOP") %>%
  select(AOP, Relevant, Rationale, `Abbreviated AOP description`, `Tox Cast Endpoints`=endPoint) %>%
  arrange(AOP) %>%
  distinct()

AOP_full_info_prior <- AOP_full_info %>%
  filter(`Tox Cast Endpoints` %in% priority_endpoints )

AOP_full_info_prior$`Tox Cast Endpoints`[sapply(AOP_full_info_prior$`Tox Cast Endpoints`, length) == 0] <- ""

y <- AOP_full_info_prior %>%
  group_by(AOP, Relevant, Rationale, `Abbreviated AOP description`) %>%
  summarise(`Tox Cast Endpoints` = list(`Tox Cast Endpoints`[`Tox Cast Endpoints` %in% priority_endpoints])) %>%
  filter(!is.na(Relevant)) 

y$`Tox Cast Endpoints`[sapply(y$`Tox Cast Endpoints`, length) == 0] <- ""

y <- y[!duplicated(y$AOP),]

fwrite(y, file = "tables/SI6.csv", na = "")
rm(list=ls())

#########################################
## SI: 9:
source(file = "data_setup.R")
AOP_crosswalk <- fread("AOP_crosswalk_Dec_2018.csv")
AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")
ear_thresh <- 0.001
siteThres <- 10

chemicalSummary <- chemicalSummary %>%
  left_join(select(AOP_crosswalk, 
                   endPoint=`Component Endpoint Name`, AOP = `AOP #`), by="endPoint") %>%
  left_join(select(AOP_info, AOP , AOP_Class = X__1), by=c("AOP"))

chemicalSummary$AOP[is.na(chemicalSummary$AOP)] <- "None"
chemicalSummary$AOP_Class[is.na(chemicalSummary$AOP_Class)] <- "Not defined"
chemicalSummary$AOP <- factor(chemicalSummary$AOP)

end_cols <- c("Not defined","Not environmentally relevant")
chemicalSummary$AOP_Class <- factor(chemicalSummary$AOP_Class, levels = c(unique(chemicalSummary$AOP_Class)[!(unique(chemicalSummary$AOP_Class) %in% end_cols)], end_cols))

chem_sum_table <- chemicalSummary %>%
  group_by(endPoint, AOP, AOP_Class, chnm,  site, date) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(endPoint, AOP, AOP_Class, chnm,  site) %>%
  summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  group_by(endPoint, AOP, AOP_Class, chnm) %>%
  summarize(hits = sum(maxEAR > ear_thresh)) %>%
  filter(hits>0) %>%
  ungroup()

chem_sum_wide <- chem_sum_table %>%
  spread(chnm, hits)

tableData2 <- select(chem_sum_wide, -endPoint, -AOP, -AOP_Class)
chem_sum_wide$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))

chem_sum_wide <- chem_sum_wide[,c("AOP", "AOP_Class", "endPoint", "nChems",
                          rev(levels(chemicalSummary$chnm))[rev(levels(chemicalSummary$chnm)) %in% names(chem_sum_wide)])]

chem_sum_wide <- chem_sum_wide %>%
  distinct() %>%
  arrange(AOP_Class, AOP)

write.csv(chem_sum_wide, file = "tables/SI9_counts.csv", row.names = FALSE, na = "")
rm(list=ls())

#########################################
# New SI 8:
source(file = "Table_SI8_fig4_Data_summaries_for_manuscript.R")
#
# source(file = "data_setup.R")
# AOP_crosswalk <- fread("AOP_crosswalk.csv")
# 
# AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")
# 
# relevance <- fread("AOP_relevance.csv") %>%
#   select(-`Endpoint(s)`) %>%
#   distinct
# 
# ear_thresh <- 0.001
# siteThres <- 10
# 
# chemicalSummary <- chemicalSummary %>%
#   left_join(select(AOP_info, endPoint = `Endpoint(s)` , AOP , AOP_Class = X__1), by=c("endPoint"))
# 
# chemicalSummary$AOP[is.na(chemicalSummary$AOP)] <- "None"
# chemicalSummary$AOP_Class[is.na(chemicalSummary$AOP_Class)] <- "Not defined"
# chemicalSummary$AOP <- factor(chemicalSummary$AOP)
# 
# end_cols <- c("Not defined","Not environmentally relevant")
# chemicalSummary$AOP_Class <- factor(chemicalSummary$AOP_Class, levels = c(unique(chemicalSummary$AOP_Class)[!(unique(chemicalSummary$AOP_Class) %in% end_cols)], end_cols))
# 
# chemicalSummary_mixtures <- chemicalSummary %>%
#   left_join(select(tox_list$chem_site, shortName=`Short Name`, Lake=site_grouping),by="shortName") %>%
#   group_by(AOP, AOP_Class, shortName, Lake, date, chnm) %>%
#   summarize(sumEAR = sum(EAR, na.rm = TRUE),
#             endpoint_count = length(unique(endPoint))) %>%
#   group_by(AOP, AOP_Class, shortName, Lake) %>%
#   summarize(maxEAR = max(sumEAR, na.rm = TRUE),
#             maxNumberEAR = max(endpoint_count),
#             chemicals = list(unique(as.character(chnm[sumEAR > ear_thresh])))) %>%
#   arrange(AOP_Class, AOP)
# 
# fwrite(chemicalSummary_mixtures, file ="tables/SI8.csv")
# rm(list=ls())
