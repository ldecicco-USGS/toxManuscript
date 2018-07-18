library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)

source(file = "data_setup.R")

hit_threshold <- 10^-3

tableData <- chemicalSummary %>%
  group_by(site, date, chnm) %>% 
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, chnm) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  group_by(chnm) %>%
  summarize(nSites = sum(meanEAR > hit_threshold)) %>%
  data.frame() %>%
  arrange(desc(nSites)) %>%
  filter(nSites > 1)

tableData$chnm <- factor(tableData$chnm, levels = tableData$chnm)

chemPlot <- ggplot(tableData)+
  geom_bar(aes(x=chnm, y=nSites),stat = "identity",fill = "steelblue") +
  theme_bw() +
  geom_text(aes(x=chnm, y=nSites-1, label = nSites), color = "white") +
  xlab("") +
  ylab("Number of Sites\n with EARmax > 0.001") +
  theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) 

chemPlot

dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(chemPlot, filename = "plots/SI3_chem_counts.png", width = 5, height = 5)
