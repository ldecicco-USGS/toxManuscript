library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)

source(file = "data_setup.R")
# Get benchmarks:
source(file = "data_setup_wq_benchmarks.R")

hit_threshold_tox <- 0.001
hit_threshold_wq <- 0.1

cas_key <- toxEval::tox_chemicals %>%
  select(CAS = Substance_CASRN, chnm_tox = Substance_Name) %>%
  full_join(select(tox_list$chem_info, CAS, orig_name = `Chemical Name`), by="CAS") %>%
  filter(CAS %in% tox_list$chem_info$CAS) 

cas_key$chnm_tox[is.na(cas_key$chnm_tox)] <- cas_key$orig_name[is.na(cas_key$chnm_tox)]
cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol monoethoxylate, (sum of all isomers; NP1EO)"] <- "4-Nonylphenol monoethoxylate" 
cas_key$chnm_tox[cas_key$chnm_tox =="4-Nonylphenol diethoxylate  (sum of all isomers; NP2EO)"] <- "4-Nonylphenol diethoxylate" 
cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol diethoxylate (OP2EO)"] <- "4-tert-Octylphenol diethoxylate" 
cas_key$chnm_tox[cas_key$chnm_tox =="4-tert-Octylphenol monoethoxylate (OP1EO)"] <- "4-tert-Octylphenol monoethoxylate" 
cas_key$chnm <- cas_key$chnm_tox

detections <- tox_list$chem_data %>%
  select(CAS, Value, SiteID) %>%
  group_by(CAS, SiteID) %>%
  top_n(n = 1, wt = Value) %>%
  filter(Value > 0) %>%
  group_by(CAS) %>%
  summarise(nSites = length(unique(SiteID))) 

tableData <- chemicalSummary %>%
  mutate(chnm = as.character(chnm)) %>%
  group_by(site, date, CAS) %>% 
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, CAS) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  full_join(detections, by="CAS") %>%
  group_by(CAS, nSites) %>%
  summarize(nSites_thres = sum(meanEAR > hit_threshold_tox, na.rm = TRUE),
            Detected = ifelse(nSites_thres>0, nSites - nSites_thres,nSites),
            has_info = ifelse(any(is.na(meanEAR)),"*","")) %>%
  ungroup() %>%
  arrange(desc(nSites_thres)) %>%
  rename(EAR = nSites_thres) %>%
  select(-nSites) %>%
  gather(criteria, number,  -CAS, -has_info) %>%
  mutate(guide_side = "ToxCast") %>%
  left_join(select(cas_key,CAS, chnm), by="CAS")

tableData_wq <- chemicalSummary_wqp %>%
  mutate(chnm = as.character(chnm)) %>%
  group_by(site, date, CAS) %>% 
  summarize(sumEAR = sum(EAR)) %>%
  group_by(site, CAS) %>%
  summarize(meanEAR = max(sumEAR)) %>%
  full_join(detections, by="CAS") %>%
  group_by(CAS, nSites) %>%
  summarize(nSites_thres = sum(meanEAR > hit_threshold_wq, na.rm = TRUE),
            Detected = ifelse(nSites_thres>0, nSites - nSites_thres,nSites),
            has_info = ifelse(any(is.na(meanEAR)),"*","")) %>%
  ungroup() %>%
  arrange(desc(nSites_thres)) %>%
  rename(TQ = nSites_thres) %>%
  select(-nSites) %>%
  gather(criteria, number,  -CAS, -has_info) %>%
  mutate(guide_side = "Water Quality Benchmarks") %>%
  left_join(select(cas_key,CAS, chnm), by="CAS")

all_data <- bind_rows(tableData, tableData_wq) 

chems_to_show <- all_data %>%
  filter(criteria %in% c("EAR","TQ"),
         number > 0) %>%
  pull(CAS)

all_data <- all_data %>%
  filter(CAS %in% chems_to_show)

all_data$chnm <- factor(all_data$chnm, levels = c(rev(unique(tableData_wq$chnm[!(unique(tableData_wq$chnm) %in% unique(tableData$chnm))])),
                                              rev(unique(tableData$chnm))))

textData <- data.frame(guide_side = c("ToxCast","Water Quality Benchmarks"),
                       textExplain = c("A","B"),
                       y = c(45,45),
                       chnm = factor(c("Metolachlor","Metolachlor"),levels = levels(all_data$chnm)),
                       stringsAsFactors = FALSE)


chemPlot <- ggplot(data = all_data) +
  geom_bar(aes(x=chnm, y=number, fill = criteria),stat = "identity") +
  geom_text(data = textData, 
            aes(x=chnm, y=y, label = textExplain)) +
  geom_point(data = filter(all_data, has_info == "*"),
             aes(x=chnm, y= 2, shape = has_info)) +
  theme_bw() +
  facet_grid(. ~ guide_side) +
  xlab("") +
  scale_fill_manual(labels = c("Detected", 
                               bquote(EAR[max] > 10^-3),
                               bquote(TQ[max] > 0.1)), 
                    values = c("#E2E2E2","#023FA5","#8E063B")) +
  scale_shape_manual(labels = "No benchmark", values = 19) +
  ylab("Number of Sites") +
  coord_flip() +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.15),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.75, "lines"),
        legend.background = element_rect(color = "black", size = 0.5, linetype = "solid")) 

chemPlot_w_cap <- chemPlot  +
  labs(caption = bquote(atop(bold("Figure SI-3:")~"Number of sites with at least one sample that resulted in an exposure activity ratio"~(EAR[SiteChem])~
                        ">"~10^-3~ 
                        "(A)","or toxicity quotient" ~ (TQ[max]) ~ "> 0.1 (B) for chemicals measured in water samples at Great Lakes tributaries, 2010-2013."))) +
  theme(plot.caption = element_text(hjust = 0, size = 7))


dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(chemPlot_w_cap, filename = "plots/SI3_chem_counts_2.pdf", width = 8, height = 5)
ggsave(chemPlot_w_cap, filename = "plots/SI3_chem_counts_2.png", width = 8, height = 5)
