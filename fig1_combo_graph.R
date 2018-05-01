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

graphData_wq <- graph_chem_data(chemicalSummary_wqp, mean_logic = "noSum")
graphData_wq$guide_side <- "Traditional\nMaximum Toxicity Quotient per Site"

graphData_eeq <- graph_chem_data(chemicalSummary_eeq)
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

dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(toxPlot_wq, filename = "plots/fig1.png", width = 12, height = 12)

