

# Figure concept for SI:
#   
# 1. Heat map with sites on the x-axis, organized like 
# the bar chart tab. AOP on the y-axis. color the heat 
# map by EAR like is done in Fig 4. Limit to AOPs that 
# have EAR > 10^-3 for at least 10 sites to keep it focused 
# on the priorities.
# 
# 2. Heat map with sites on the x-axis, organized like 
# the bar chart tab. Chemicals on the y-axis. Color the heat 
# map by EAR like is done in Fig 4. Limit to chemicals 
# that have EAR > 10^-3 to keep it sane.
# --start this one using the Heat map tab from the app, but 
# reducing the chemicals that have EAR > 10^-3 to keep it focused. 
# We may or may not need the categories on the right, but leave 
# them in for now.

library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)
####################################
source(file = "data_setup.R")

ear_thresh <- 0.001
siteThres <- 10

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

boxData_max <- chemicalSummary %>%
  left_join(AOP, by="endPoint") %>%
  group_by(ID, chnm, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)])

boxData <- boxData_max %>%
  group_by(ID,site,date) %>%
  summarize(total = sum(maxEAR))  %>%
  group_by(ID, site) %>%
  summarize(maxMaxEAR = max(total, na.rm = TRUE),
            date_picked = date[which.max(total)]) %>%
  ungroup() %>%
  filter(!is.na(ID)) 

priority_AOPs <- boxData %>%
  filter(maxMaxEAR > ear_thresh) %>%
  group_by(ID) %>%
  summarise(siteDet = n_distinct(site)) %>%
  filter(siteDet >= siteThres)

boxData_max <- filter(boxData_max, ID %in% priority_AOPs$ID)
boxData <- filter(boxData, ID %in% priority_AOPs$ID)

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)

relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.) %>%
  filter(ID %in% priority_AOPs$ID) 

boxData <- boxData %>%
  left_join(select(relevance, ID, Relevant), by="ID") %>%
  mutate(ID = factor(ID)) 

boxData <- boxData %>%
  left_join(select(tox_list$chem_site, site=SiteID, site_grouping, name=`Short Name`), by="site")

boxData <- filter(boxData, Relevant %in% c("Yes","Maybe"))

si1 <- ggplot(data = boxData) +
  geom_tile(aes(x = name, y = ID, fill = maxMaxEAR)) +
  theme_bw() +
  ylab("AOP ID") +
  labs(fill="Max EAR") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid( ~ site_grouping, scales="free", space="free") +
  scale_y_discrete(drop=TRUE) +
  scale_fill_gradient( guide = "legend",
                       trans = 'log', limits = c(1e-8,1),
                       low = "white", high = "steelblue",
                       breaks = c(1e-8,1e-6,1e-4,1e-2,1),
                       labels = toxEval:::fancyNumbers2,
                       na.value = 'transparent') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.ticks = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill="transparent"),
        panel.spacing = unit(0.05, "lines"),
        plot.background = element_rect(fill = "transparent",colour = NA))

png("plots/si3.png", width = 1200, height = 800, res = 142)
si1
dev.off()


#################################################3
# Now by chemical:
source(file = "plot_tox_heatmap_manuscript.R")


plot_out <- plot_tox_heatmap_manuscript(chemicalSummary, tox_list$chem_site,  category = "Chemical")

ggsave(plot_out, file = "plots/si4.png", width = 11, height = 9)

