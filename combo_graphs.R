library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################
source(file = "data_setup.R")

####################################
# Get benchmarks:
source(file = "data_setup_wq_benchmarks.R")

####################################
# Get concentrations:
source(file = "data_setup_concentrations.R")

# Special funtion:
source(file = "combo_graph_function.R")

graphData_tox <- graph_chem_data(chemicalSummary)
graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"

toxPlot_wq <- combo_plot_matches(graphData_tox, graphData_wq, thres_1 = 10^-3, thres_2 = 10^-1, drop = FALSE)

gb <- ggplot2::ggplot_build(toxPlot_wq)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-2-1"] <- "off"

png("WQ_Tox_no_drop.png", width = 1200, height = 1200, res = 142)
grid::grid.draw(gt)
dev.off()

toxPlot_wq_no_thres <- combo_plot_matches(graphData_tox, graphData_wq, thres_1 = NA, thres_2 = NA, drop = FALSE)

gb <- ggplot2::ggplot_build(toxPlot_wq_no_thres)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-2-1"] <- "off"

png("WQ_Tox_no_drop_no_thres.png", width = 1200, height = 1200, res = 142)
grid::grid.draw(gt)
dev.off()

toxPlot_wq_drop <- combo_plot_matches(graphData_tox, graphData_wq, thres_1 = 10^-3, thres_2 = 10^-1, drop = TRUE)

gb <- ggplot2::ggplot_build(toxPlot_wq_drop)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-2-1"] <- "off"

png("WQ_Tox.png", width = 1200, height = 900, res = 142)
grid::grid.draw(gt)
dev.off()

toxPlot_conc <- combo_plot_matches(graphData_tox, graphData_conc, thres_1 = NA, thres_2 = NA, drop = FALSE)

gb <- ggplot2::ggplot_build(toxPlot_conc)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-2-1"] <- "off"

png("Conc_Tox_no_drop.png", width = 1200, height = 1200, res = 142)
grid::grid.draw(gt)
dev.off()

toxPlot_eeq <- combo_plot_matches(graphData_tox, graphData_eeq, thres_1 = NA, thres_2 = NA)

gb <- ggplot2::ggplot_build(toxPlot_eeq)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"
gt$layout$clip[gt$layout$name=="panel-2-1"] <- "off"

dir.create(file.path("plots"), showWarnings = FALSE)
png("plots/EEQ_Toxp.png", width = 1000, height = 400, res = 142)
grid::grid.draw(gt)
dev.off()

