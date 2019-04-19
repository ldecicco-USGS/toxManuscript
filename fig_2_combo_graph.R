library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

####################################
# Get benchmarks:
source(file = "data_setup_wq_benchmarks.R")

# Special funtion:
source(file = "combo_graph_function.R")

# guide_side is the title of the side-by-side labels
# or...the column headers
graphData_tox <- graph_chem_data_CAS(chemicalSummary)
graphData_tox$guide_side <- "atop(ToxCast,Maximum~EAR[SiteChem])"

graphData_wq <- graph_chem_data_CAS(chemicalSummary_wqp, sum_logic = FALSE)
graphData_wq$guide_side <- "atop(Traditional,Maximum~Toxicity~Quotient~per~Site)"

graphData_eeq <- graph_chem_data_CAS(chemicalSummary_eeq, sum_logic = FALSE)
graphData_eeq$guide_side <- "atop(Traditional,Maximum~Toxicity~Quotient~per~Site)"

# guide_up is the top-and-bottom labels. We're currently not showing those:

graphData_wq$guide_up <- "A"
graphData_eeq$guide_up <- "B"
graphData_tox$guide_up <- "A"

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

# cas_key$chnm <- paste0(cas_key$chnm_tox, "\t(",cas_key$nSites,")")

# Duplicating the "main" data let's us use the "drop" and astrict calcs:
graphData_tox_dup <- graphData_tox %>%
  filter(CAS %in% graphData_eeq$CAS) %>%
  mutate(guide_up = "B")

graphData_tox <- bind_rows(graphData_tox, graphData_tox_dup)

graphData_tox <- graphData_tox %>%
  left_join(select(cas_key, CAS, chnm), by="CAS")

graphData_wq <- graphData_wq %>%
  left_join(select(cas_key, CAS, chnm), by="CAS")

graphData_eeq <- graphData_eeq %>%
  left_join(select(cas_key, CAS, chnm), by="CAS")

# drop is whether or not to drop empty rows
# grid is if we want it in a 2x2 grid (a la WQB and EEQ) or
# (false) 3 rows (a la tox/WQB/Conc)
toxPlot_wq <- combo_plot_matches(graphData_tox, graphData_wq, 
                                 thres_1 = NA, thres_2 = NA, 
                                 drop = FALSE, grid = TRUE, 
                                 gd_3 = graphData_eeq)

plot_data <- toxPlot_wq$data

textData <- data.frame(guide_up = c("A","A","B","B"),
                       guide_side = factor(c(levels(plot_data$guide_side)[1],
                                             levels(plot_data$guide_side)[2],
                                             levels(plot_data$guide_side)[1],
                                             levels(plot_data$guide_side)[2]),levels = levels(plot_data$guide_side))) %>%
  mutate(textExplain = c("A","B","C","D"),
         y = c(10,100,10,100),
         chnm = factor(c("4-Nonylphenol, branched","4-Nonylphenol, branched",
                         "4-tert-Octylphenol monoethoxylate","4-tert-Octylphenol monoethoxylate"), 
                       levels = levels(plot_data$chnm)))

toxPlot_wq <- toxPlot_wq +
  geom_text(data = textData, aes(label = textExplain, x = chnm, y=y),size = 3)


no_axis <- toxPlot_wq +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

site_counts <- tox_list$chem_data %>%
  filter(Value > 0) %>%
  group_by(CAS) %>%
  summarise(nSites = length(unique(SiteID))) %>%
  ungroup() %>%
  right_join(distinct(select(plot_data, CAS, Class, chnm, guide_up)), by="CAS") %>%
  mutate(guide_side = "atop(N,Sites[ ])")

site_counts$nSites[is.na(site_counts$nSites)] <- 0

pretty_logs_new <- toxEval:::prettyLogs(c(10,100,1000))

site_graph <- ggplot(data = site_counts) +
  geom_text(aes(x=chnm, label = nSites, y=100), size = 3) +
  facet_grid(guide_up ~ guide_side, scales = "free", space = "free_y", labeller = label_parsed)+
  theme_bw() +
  coord_flip() +
  scale_y_log10(labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
  theme(axis.text.x = element_text( color = "transparent"),
        axis.text.y = element_text(size=9, vjust=.35, color = "black"),
        axis.title=element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill = "transparent",colour = NA),
        strip.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "transparent")) 

library(cowplot)
legend_box <- get_legend(no_axis + 
                           theme(legend.position = "bottom") )

dir.create(file.path("plots"), showWarnings = FALSE)

png("plots/Fig2.png", width = 1800, height = 1200, res = 142)
plot_grid(site_graph, no_axis,
          NULL,legend_box,
          align = "h", nrow = 2,ncol=2, 
          rel_widths =  c(2/9, 7/9),
          rel_heights = c(0.9,0.1))

dev.off()

pdf("plots/Fig2.pdf", width = 11, height = 9)
plot_grid(site_graph, no_axis,
          NULL,legend_box,
          align = "h", nrow = 2,ncol=2, 
          rel_widths =  c(2.5/9, 6.5/9),
          rel_heights = c(0.9,0.1))

dev.off()

# ggsave(toxPlot_wq, filename = "plots/fig1.png", width = 12, height = 12)

# dir.create("plots", showWarnings = FALSE)
# png("plots/fig2_no_clip_new.png", width = 1200, height = 1200, res = 142)
# gb <- ggplot2::ggplot_build(toxPlot_wq)
# gt <- ggplot2::ggplot_gtable(gb)
# 
# gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
# gt$layout$clip[gt$layout$name=="panel-1-2"] <- "off"
# grid::grid.draw(gt)
# dev.off()
