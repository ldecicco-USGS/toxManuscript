library(dplyr)
library(tidyr)
library(readxl)

source(file = "data_setup.R")

ear_thresh <- 0.001
siteThres <- 10

## SI 1&2

# si_1 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-1 Watershed Characteristics")
# si_2 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-2 Chemical Classes")

# SI 3:
# 
# End Point|End point long name|Assay Source| AOP ID | Relavence

endPointInfo <- clean_endPoint_info(endPointInfo)

si_3_endpoints <- select(endPointInfo, 
                         `ToxCast Endpoint`=assay_component_endpoint_name, 
                         `Assay Source`=assay_source_long_name) %>%
  distinct() 

dir.create("tables", showWarnings = FALSE)
write.csv(si_3_endpoints, file = "tables/SI3.csv", row.names = FALSE, na = "")

# SI 4:

#look for code already in sup, remove redundent stuff

# SI 5:
#
# 
library(data.table)
AOP_crosswalk <- fread("AOP_crosswalk.csv")
write.csv(AOP_crosswalk, file = "tables/SI5.csv", row.names = FALSE, na = "")


# SI 6:
#
# 
relevance <- fread("AOP_relevance.csv") %>%
  select(-`Endpoint(s)`) %>%
  distinct
write.csv(relevance, file = "tables/SI6.csv", row.names = FALSE, na = "")

# chemicalSummary <- chemicalSummary %>%
#   left_join(select(endPointInfo, 
#                    endPoint=assay_component_endpoint_name,
#                    subFamily = intended_target_family_sub), by="endPoint") %>%
#   left_join(AOP, by="endPoint")
# 
# nSite_table <- chemicalSummary %>%
#   rename(Family = Bio_category) %>%
#   group_by(endPoint, chnm, subFamily, Family, ID, site, date) %>%
#   summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
#   group_by(endPoint, chnm, subFamily, Family, ID, site) %>%
#   summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
#   group_by(endPoint, chnm, subFamily, Family, ID) %>%
#   filter(n_distinct(site) > siteThres) %>%
#   summarize(hits = sum(maxEAR > ear_thresh)) %>%
#   filter(hits>0)
# 
# nSite_table_wide <- nSite_table %>%
#   spread(chnm, hits) %>%
#   arrange(Family, subFamily, ID, endPoint) %>%
#   select(Family, subFamily, ID, endPoint, everything()) %>%
#   ungroup()
# 
# tableData2 <- select(nSite_table_wide, -endPoint, -Family, -subFamily, -ID)
# nSite_table_wide$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))
# 
# nSite_table_wide <- nSite_table_wide[,c("Family", "subFamily","ID",
#                           "endPoint", "nChems",
#                           rev(levels(chemicalSummary$chnm))[rev(levels(chemicalSummary$chnm)) %in% names(nSite_table_wide)])]  
# 
# dir.create("tables", showWarnings = FALSE)
# write.csv(nSite_table_wide, file = "tables/SI5.csv", row.names = FALSE, na = "")

# New SI 6:
#
# Site | AOP ID | AOP name | EndPoint | Chemical | max EAR | mean EAR
# rm(chemicalSummary)
# source(file = "data_setup.R")
# 
# si_6 <- chemicalSummary %>%
#   left_join(select(AOP_crosswalk, 
#                    endPoint=Component.Endpoint.Name, 
#                    ID=AOP..,
#                    AOP_name = AOP.Title), by="endPoint") %>%
#   left_join(select(tox_list$chem_site, site=SiteID, `Short Name`),by="site") %>%
#   group_by(ID, AOP_name, `Short Name`, date, chnm, endPoint) %>%
#   summarize(sumEAR = sum(EAR, na.rm = TRUE)) %>%
#   group_by(ID, AOP_name, `Short Name`, chnm, endPoint) %>%
#   summarise(meanEAR = mean(sumEAR, na.rm = TRUE),
#             maxEAR = max(sumEAR, na.rm = TRUE)) %>%
#   ungroup() %>%
#   filter(meanEAR > ear_thresh) %>%
#   rename(site=`Short Name`)
# 
# write.csv(si_6, file = "tables/SI6.csv", row.names = FALSE, na = "")
# 
# 
# # SI 6:
# #
# # 
