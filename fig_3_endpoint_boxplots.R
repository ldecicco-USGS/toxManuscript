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
siteThreshold <- 20

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= threshold) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites))

siteRows <- which(endpoints_sites_hits$numSites >= siteThreshold)
priority_endpoints <- pull(endpoints_sites_hits,endPoint)[siteRows]

chemicalSummaryPriority <- chemicalSummary[which(chemicalSummary$endPoint %in% priority_endpoints),]

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

eps_with_ids <- unique(AOP$endPoint)

chemicalSummaryPriority$has_AOP <- "AOP Undefined"
chemicalSummaryPriority$has_AOP[chemicalSummaryPriority$endPoint %in% eps_with_ids] <- "AOP Associated"

endpointPlot <- plot_tox_endpoints_manuscript(chemicalSummaryPriority)

#  not sure how to put the number of chemicals per endpoint on the right side.
#  need to add color to distinguish endpoints with associated AOPs

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

sitesChemsPerEndoint <- left_join(endpoints_unique_chems,endpoints_sites_hits)

# unique(test$CAS)
# ATG_PXRE_CIS_up
# CLD_CYP1A1_6hr
# write.csv(endpoints_sites_hits,file="sitesChemsPerEndoint.csv",row.names = FALSE)
