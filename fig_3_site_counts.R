library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

threshold <- 10^-3

graphData <- chemicalSummary %>%
  group_by(site, date, chnm) %>% 
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, chnm) %>%
  summarize(maxEAR = max(sumEAR),
            count = n()) %>%
  group_by(site, count) %>%
  summarize(nChem = sum(maxEAR > threshold)) %>%
  data.frame() %>%
  left_join(select(tox_list[["chem_site"]], site=SiteID, `Short Name`, site_grouping),
            by = "site") %>%
  mutate(`Short Name` = factor(`Short Name`, levels = sitesOrdered)) %>%
  mutate(site_grouping = factor(site_grouping, levels = c("Lake Superior",
                                                          "Lake Michigan",
                                                          "Lake Huron",
                                                          "Lake Erie",
                                                          "Lake Ontario")))

placement <- -0.05*diff(range(graphData$nChem))
chem_site <- tox_list[["chem_site"]]
label_samples <- data.frame(x=-Inf,
                            y=placement,
                            label="# Samples", 
                            site_grouping = NA,
                            stringsAsFactors = FALSE)
label_samples$site_grouping <- factor(levels(chem_site$site_grouping)[1],
                                      levels = levels(chem_site$site_grouping))

countPlot <- ggplot(graphData, aes(x=`Short Name`))+
  geom_bar(aes(y=nChem),
           stat = "identity",
           fill = "steelblue") +
  geom_text(aes(y=-1, label =  count), size = 2.5, angle = 90) +
  theme_bw() +
  facet_grid(. ~ site_grouping, scales="free", space="free") +
  xlab("") +
  ylab("Number of Chemicals") +
  theme(strip.text.y = element_text(angle=0, hjust=0, size=7), 
        strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill="transparent", colour = NA),
        axis.text = element_text(size=8),
        axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1),
        panel.spacing = unit(0.15, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  geom_text(data = label_samples,vjust=0.5,hjust=1.1,
            aes(x=x,y=y,label=label),
            size=2,inherit.aes = FALSE) +
  coord_cartesian(clip="off")

countPlot



dir.create(file.path("plots"), showWarnings = FALSE)
png("plots/fig3_site_count.png", width = 1000, height = 800, res = 142)
print(countPlot)
dev.off()
