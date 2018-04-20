source(file = "data_setup.R")

stack_plot <- plot_tox_stacks(chemicalSummary, 
                chem_site = tox_list$chem_site,
                category = 'Biological',
                mean_logic = FALSE,
                include_legend = FALSE)

gb <- ggplot2::ggplot_build(stack_plot)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel-1-1"] <- "off"

png("stacks.png", width = 1100, height = 700, res = 142)
grid::grid.draw(gt)
dev.off()

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                              category = "Biological") 
cbValues <- colorRampPalette(cbPalette)(length(levels(graphData$category)))
set.seed(4)
cbValues <- sample(cbValues)
names(cbValues) <- levels(graphData$category)

bp_plot <- plot_tox_boxplots(chemicalSummary,
                             category = "Biological", 
                             pallette = cbValues, 
                             font_size = 20)

gb <- ggplot2::ggplot_build(bp_plot)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel"] <- "off"

png("box_bio.png", width = 1000, height = 700, res = 142)
grid::grid.draw(gt)
dev.off()

