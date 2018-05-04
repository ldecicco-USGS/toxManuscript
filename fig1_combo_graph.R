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
graphData_tox <- graph_chem_data(chemicalSummary)
graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"

graphData_wq <- graph_chem_data(chemicalSummary_wqp, sum_logic = FALSE)
graphData_wq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"

graphData_eeq <- graph_chem_data(chemicalSummary_eeq, sum_logic = FALSE)
graphData_eeq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"

# guide_up is the top-and-bottom labels. We're currently not showing those:

graphData_wq$guide_up <- "A"
graphData_eeq$guide_up <- "B"
graphData_tox$guide_up <- "A"

# Duplicating the "main" data let's us use the "drop" and astrict calcs:
graphData_tox_dup <- graphData_tox
graphData_tox_dup$guide_up <- "B"

graphData_tox_dup <- filter(graphData_tox_dup, graphData_tox_dup$chnm %in% graphData_eeq$chnm)
graphData_tox <- bind_rows(graphData_tox, graphData_tox_dup)

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
                         "4-tert-Octylphenol monoethoxylate (OP1EO)","4-tert-Octylphenol monoethoxylate (OP1EO)"), 
                       levels = levels(plot_data$chnm)))

toxPlot_wq <- toxPlot_wq +
  geom_text(data = textData, aes(label = textExplain, x = chnm, y=y),size = 3)

dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(toxPlot_wq, filename = "plots/fig1.png", width = 12, height = 12)


png("plots/fig1_no_clip.png", width = 1200, height = 1200, res = 142)
gb <- ggplot2::ggplot_build(toxPlot_wq)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-1-2"] <- "off"
grid::grid.draw(gt)
dev.off()
