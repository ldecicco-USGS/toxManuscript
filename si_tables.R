library(dplyr)
library(tidyr)

source(file = "data_setup.R")

ear_thresh <- 0.001
siteThres <- 10

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct() 

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)
relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.)  

# SI 3:
# 
# End Point|End point long name|Assay Source| AOP ID | Relavence

endPointInfo <- clean_endPoint_info(endPointInfo)

si_3_endpoints <- AOP %>%
  left_join(select(relevance, ID, Relevant), by="ID") %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name, 
                   assay_source=assay_source_long_name),
            by="endPoint") 
dir.create("tables", showWarnings = FALSE)
write.csv(si_3_endpoints, file = "tables/SI3.csv", row.names = FALSE, na = "")

# SI 5:
#
# 

chemicalSummary <- chemicalSummary %>%
  left_join(select(endPointInfo, 
                   endPoint=assay_component_endpoint_name,
                   subFamily = intended_target_family_sub), by="endPoint") %>%
  left_join(AOP, by="endPoint")

nSite_table <- chemicalSummary %>%
  rename(Family = Bio_category) %>%
  group_by(endPoint, chnm, subFamily, Family, ID, site, date) %>%
  summarise(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(endPoint, chnm, subFamily, Family, ID, site) %>%
  summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  group_by(endPoint, chnm, subFamily, Family, ID) %>%
  filter(n_distinct(site) > siteThres) %>%
  summarize(hits = sum(maxEAR > ear_thresh)) %>%
  filter(hits>0)

nSite_table_wide <- nSite_table %>%
  spread(chnm, hits) %>%
  arrange(Family, subFamily, ID, endPoint) %>%
  select(Family, subFamily, ID, endPoint, everything()) %>%
  ungroup()

tableData2 <- select(nSite_table_wide, -endPoint, -Family, -subFamily, -ID)
nSite_table_wide$nChems <- apply(tableData2, MARGIN = 1, function(x) sum(x>0, na.rm = TRUE))

nSite_table_wide <- nSite_table_wide[,c("Family", "subFamily","ID",
                          "endPoint", "nChems",
                          rev(levels(chemicalSummary$chnm))[rev(levels(chemicalSummary$chnm)) %in% names(nSite_table_wide)])]  

dir.create("tables", showWarnings = FALSE)
write.csv(nSite_table_wide, file = "tables/SI5.csv", row.names = FALSE, na = "")

# New SI 6:
#
# Site | AOP ID | AOP name | EndPoint | Chemical | max EAR | mean EAR
rm(chemicalSummary)
source(file = "data_setup.R")

si_6 <- chemicalSummary %>%
  left_join(select(AOP_crosswalk, 
                   endPoint=Component.Endpoint.Name, 
                   ID=AOP..,
                   AOP_name = AOP.Title), by="endPoint") %>%
  left_join(select(tox_list$chem_site, site=SiteID, `Short Name`),by="site") %>%
  group_by(ID, AOP_name, `Short Name`, date, chnm, endPoint) %>%
  summarize(sumEAR = sum(EAR, na.rm = TRUE)) %>%
  group_by(ID, AOP_name, `Short Name`, chnm, endPoint) %>%
  summarise(meanEAR = mean(sumEAR, na.rm = TRUE),
            maxEAR = max(sumEAR, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(meanEAR > ear_thresh) %>%
  rename(site=`Short Name`)

write.csv(si_6, file = "tables/SI6.csv", row.names = FALSE, na = "")
