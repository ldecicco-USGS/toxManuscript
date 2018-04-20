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

graphData_tox <- graph_chem_data(chemicalSummary)
graphData_tox$guide_side <- "ToxCast\nMaximum EAR per Site"

combo_plot_matches <- function(gd_1, gd_2,
                               thres_1, thres_2,
                               drop = TRUE){

  orderChem_1_2 <- bind_rows(gd_1, 
                             filter(gd_2, !(chnm %in% levels(gd_1$chnm)))) %>%
    group_by(chnm,Class) %>%
    summarise(median = quantile(maxEAR[maxEAR != 0],0.5)) %>%
    data.frame() 
  
  class_order <- toxEval:::orderClass(bind_rows(gd_1, 
                                                filter(gd_2, !(chnm %in% levels(gd_1$chnm)))))
  
  orderChem_1_2 <- orderChem_1_2 %>%
    mutate(Class = factor(Class, levels=class_order$Class)) %>%
    arrange(desc(Class), desc(!is.na(median)), median)
  
  graphData_1_2 <- bind_rows(gd_1, gd_2)
  graphData_1_2$Class <- factor(graphData_1_2$Class, levels = class_order$Class)
  graphData_1_2$chnm <- factor(graphData_1_2$chnm, levels = orderChem_1_2$chnm)
  
  chems_in_1_2 <- select(graphData_1_2, guide_side, chnm) %>%
    distinct() %>%
    select(-guide_side) %>%
    table() %>%
    data.frame() %>%
    filter(Freq == 2) %>%
    mutate(chnm = unlist(`.`)) %>%
    select(chnm)
  
  if(drop){
    graphData_1_2 <- filter(graphData_1_2, chnm %in% chems_in_1_2$chnm)
  }
  
  guide_side_1 <- gd_1$guide_side[1]
  guide_side_2 <- gd_2$guide_side[1]
  
  countNonZero_1_2 <- graphData_1_2 %>%
    group_by(site, chnm, Class, guide_side) %>%
    summarise(meanEAR = mean(maxEAR, na.rm=TRUE)) %>%
    group_by(chnm, Class, guide_side) %>%
    summarise(nonZero = as.character(sum(meanEAR>0)),
              hits = as.character(sum(meanEAR > ifelse(guide_side == guide_side_1, thres_1, thres_2)))) %>%
    data.frame() 
  
  thresh_df <- data.frame(guide_side = c(guide_side_1,guide_side_2),
                          thres = c(thres_1, thres_2),
                          stringsAsFactors = FALSE)
  
  thresh_df <- thresh_df[!is.na(thresh_df$thres),]
  
  graphData_1_2$guide_side <- factor(graphData_1_2$guide_side, levels = c(guide_side_1,guide_side_2))
  thresh_df$guide_side <- factor(thresh_df$guide_side, levels = c(guide_side_1,guide_side_2))
  countNonZero_1_2$guide_side <- factor(countNonZero_1_2$guide_side, levels = c(guide_side_1,guide_side_2))
  
  cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
  toxPlot_1_2 <- ggplot(data=graphData_1_2) +
    scale_y_log10(labels=toxEval:::fancyNumbers)  +
    geom_boxplot(aes(x=chnm, y=maxEAR, fill=Class),
                 lwd=0.1,outlier.size=1) +
    facet_grid(. ~ guide_side, scales = "free") +
    theme_minimal() +
    coord_flip() +
    theme(axis.text = element_text( color = "black"),
          axis.text.y = element_text(size=7),
          axis.title=element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank()) +
    guides(fill=guide_legend(ncol=6)) +
    theme(legend.position="bottom",
          legend.justification = "left",
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.title=element_blank(),
          legend.text = element_text(size=8),
          legend.key.height = unit(1,"line"),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = cbValues, drop=FALSE)
  
  if(!drop){
    toxPlot_1_2 <- toxPlot_1_2 +
      scale_x_discrete(drop=TRUE) 
  }
  
  plot_info <- ggplot_build(toxPlot_1_2)
  layout_stuff <- plot_info$layout
  
  xmin_1 <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
  xmax_1 <- 10^(layout_stuff$panel_ranges[[1]]$x.range[2])
  xmin_2 <- 10^(layout_stuff$panel_ranges[[2]]$x.range[1])
  xmax_2 <- 10^(layout_stuff$panel_ranges[[2]]$x.range[2])

  ymax <- floor(layout_stuff$panel_ranges[[2]]$y.range[2])

  if(!drop){
    chm_side <- bind_rows(select(orderChem_1_2,-median),select(orderChem_1_2,-median) ) %>%
      mutate(guide_side = factor(c(rep(guide_side_1,nrow(orderChem_1_2)),
                                   rep(guide_side_2,nrow(orderChem_1_2))),
                                 levels = c(guide_side_1, guide_side_2)))
    chm_side$chnm <- factor(chm_side$chnm, levels = levels(graphData_1_2$chnm))
    
    countNonZero_1_2 <- right_join(countNonZero_1_2, chm_side, by=c("chnm","Class","guide_side"))
    countNonZero_1_2$nonZero[is.na(countNonZero_1_2$nonZero)] <- "*"
  }
    
  countNonZero_1_2 <- countNonZero_1_2 %>%
    mutate(ymin = ifelse(guide_side == guide_side_1, xmin_1, xmin_2),
           ymax = ifelse(guide_side == guide_side_1, xmax_1, xmax_2))
  


  labels_1_2 <- data.frame(x = rep(Inf,4),
                           y = c(xmin_1,thres_1,
                                 xmin_2,thres_2),
                           label = rep(c("# Sites","Threshold"),2),
                           guide_side = c(rep(guide_side_1,2),rep(guide_side_2,2)),
                           stringsAsFactors = FALSE)
  
  if(!all(is.na(countNonZero_1_2$hits[countNonZero_1_2$guide_side == guide_side_1]))){
    labels_1_2 <- bind_rows(labels_1_2,
                            data.frame(x = Inf, y = xmax_1, 
                                       label = "# Hits", 
                                       guide_side = guide_side_1,
                                       stringsAsFactors = FALSE))
  }
  
  if(!all(is.na(countNonZero_1_2$hits[countNonZero_1_2$guide_side == guide_side_2]))){
    labels_1_2 <- bind_rows(labels_1_2,
                            data.frame(x = Inf, y = xmax_2, 
                                       label = "# Hits", 
                                       guide_side = guide_side_2,
                                       stringsAsFactors = FALSE))
  }
  
  labels_1_2$guide_side <- factor(labels_1_2$guide_side, levels = c(guide_side_1, guide_side_2))
  
  toxPlot_1_2 <- toxPlot_1_2 +
    geom_text(data=countNonZero_1_2, size=2.5, 
              aes(x= chnm, label = nonZero, y=ymin)) +
    geom_text(data=countNonZero_1_2, size=2.5, 
              aes(x= chnm, label = hits, y=ymax)) +
    geom_text(data=labels_1_2, size=2.5,
              aes(x = x,  y=y, label = label)) +
    geom_segment(data = thresh_df, aes(y = thres, yend = thres),
                 linetype="dashed", 
                 x = 1,
                 xend = ymax) 
  
  return(toxPlot_1_2)
  
}


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

