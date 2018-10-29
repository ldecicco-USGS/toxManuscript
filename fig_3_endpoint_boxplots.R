library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")
source(file = "plot_tox_endpoints_manuscript.R")

#1. Determine which endpoints have EARmax > threshold and at how many sites
#2. For endpoints from #1, how many and which chemicals contribute for each endpoint
#3. Develop boxplots of EARmax by endpoint with number of sites on the left and number of chems on the right
#4. Match with AOPs and color the boxplots differently for those with AOPs and those without

threshold <- 0.001
siteThreshold <- 10

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= threshold) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThreshold)

priority_endpoints <- endpoints_sites_hits$endPoint

chemicalSummaryPriority <- filter(chemicalSummary, endPoint %in% priority_endpoints)

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

eps_with_ids <- unique(AOP$endPoint)

chemicalSummaryPriority$has_AOP <- "AOP Undefined"
chemicalSummaryPriority$has_AOP[chemicalSummaryPriority$endPoint %in% eps_with_ids] <- "AOP Associated"

endpointPlot <- plot_tox_endpoints_manuscript(chemicalSummaryPriority, 
                                              category = "Chemical", 
                                              font_size = 7,title = " ",
                                              pallette = c("steelblue", "white"))

gb <- ggplot2::ggplot_build(endpointPlot)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel"] <- "off"

dir.create(file.path("plots"), showWarnings = FALSE)
png("plots/fig3_endpoint_boxplots.png", width = 1000, height = 800, res = 142)
grid::grid.draw(gt)
dev.off()

# Determine number of chemicals and sites per endpoint
endpoints_unique_chems <- filter(chemicalSummaryPriority,EAR > 0) %>%
  group_by(endPoint) %>%
  summarize(numChems = n_distinct(CAS))

endpoints_unique_sites <- filter(chemicalSummaryPriority,EAR > 0) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site))

sitesChemsPerEndoint <- left_join(endpoints_unique_chems,endpoints_unique_sites)

####Determine a few things for the text: 
#how many endpoints when using threshold and siteThreshold
# unique(chemicalSummaryPriority$endPoint)
# unique(endpoints_sites_hits$endPoint) #48 endpoints for threshold = 0.001 and siteThreshold = 10
# 
# # merge the priority endpoints with the endpoint crosswalk for transmittal to EPA for relevance evaluation
# AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
# AOP_OWC <-  left_join(AOP_crosswalk,endpoints_sites_hits,by=c("Component.Endpoint.Name"="endPoint")) %>%
#   filter(!is.na(numSites))
# unique(AOP_OWC$Assay.Endpoint.ID)
# write.csv(AOP_OWC,file="AOPs_for_Great_Lakes_OWC_study.csv",row.names = FALSE)

# How many endpoints per AOP

endpoints_per_AOP <- left_join(AOP,chemicalSummaryPriority,by="endPoint") %>%
  filter(EAR > 0) %>%
  group_by(ID) %>%
  summarize(numEndpoints = n_distinct(endPoint))
range(endpoints_per_AOP$numEndpoints)

endpoints_per_AOP <- left_join(AOP,chemicalSummary,by="endPoint") %>%
  filter(EAR > 0) %>%
  group_by(ID) %>%
  summarize(numEndpoints = n_distinct(endPoint))
range(endpoints_per_AOP$numEndpoints)
