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
  summarize(nSites_thres = sum(meanEAR > hit_threshold),
            Detected = sum(meanEAR > 0) - nSites_thres) %>%
  data.frame() %>%
  arrange(desc(nSites_thres)) %>%
  filter(nSites_thres > 1) %>%
  gather(criteria, number,  -chnm)

tableData$chnm <- factor(tableData$chnm, levels = unique(tableData$chnm))

chemPlot <- ggplot(tableData)+
  geom_bar(aes(x=chnm, y=number, fill = criteria),stat = "identity") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(labels = c("Detected", bquote(EAR[max] > 10^-3)), values = c("grey","steelblue")) +
  ylab("Number of Sites") +
  theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(color = "black", size = 0.5, linetype = "solid")) 

chemPlot_w_cap <- chemPlot  +
  labs(caption = bquote(atop(bold("Figure SI-3:")~"Number of sites with at least one sample that resulted in an exposure activity ratio"~ (EAR[SiteChem]),
                        ">" ~  10^-3 ~ 
                        "for chemicals measured in water samples at Great Lakes tributaries, 2010-2013.        "))) +
  theme(plot.caption = element_text(hjust = 0, size = 6)) 


dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(chemPlot_w_cap, filename = "plots/SI3_chem_counts.pdf", width = 5, height = 5)
