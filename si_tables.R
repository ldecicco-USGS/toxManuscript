library(dplyr)
library(tidyr)
library(readxl)
library(data.table)

source(file = "data_setup.R")

AOP_crosswalk <- fread("AOP_crosswalk.csv")

AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")

relevance <- fread("AOP_relevance.csv") %>%
  select(-`Endpoint(s)`) %>%
  distinct

ear_thresh <- 0.001
siteThres <- 10

## SI 1&2

# si_1 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-1 Watershed Characteristics")
# si_2 <- readxl::read_excel("M:/QW Monitoring Team/GLRI toxics/ToxCast/JAs/ToxEval 1/Manuscript/Supplemental.xlsx",sheet = "SI-2 Chemical Classes")

# SI 3:
# 
# End Point|End point long name|Assay Source| AOP ID | Relavence

endPointInfo <- clean_endPoint_info(end_point_info)

si_3_endpoints <- select(endPointInfo, 
                         `ToxCast Endpoint`=assay_component_endpoint_name, 
                         `Assay Source`=assay_source_long_name) %>%
  distinct() 

si_3_endpoints <- si_3_endpoints %>%
  left_join(select(AOP_info, `Endpoint(s)`, AOP, Relevant), by=c("ToxCast Endpoint"="Endpoint(s)"))

dir.create("tables", showWarnings = FALSE)
write.csv(si_3_endpoints, file = "tables/SI3.csv", row.names = FALSE, na = "")

# SI 4:

# si_table_counts.R

# SI 5:
#
# 
write.csv(AOP_crosswalk, file = "tables/SI5.csv", row.names = FALSE, na = "")


# SI 6:
#
# 
relevance <- fread("AOP_relevance.csv") %>%
  select(-`Endpoint(s)`) %>%
  distinct
write.csv(relevance, file = "tables/SI6.csv", row.names = FALSE, na = "")


## SI: 7:

# chemicalSummary <- chemicalSummary %>%
#   left_join(select(AOP_crosswalk, endPoint = `Component Endpoint Name` , `AOP ID`= `AOP #`), by=c("endPoint"))
chemicalSummary <- chemicalSummary %>%
  left_join(select(AOP_info, endPoint = `Endpoint(s)` , AOP , AOP_Class = X__1), by=c("endPoint"))

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

write.csv(chem_sum_wide, file = "tables/SI7_counts.csv", row.names = FALSE, na = "")


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
