library(toxEval)
library(dplyr)
library(tidyr)
library(stringi)

source(file = "data_setup.R")

# Substitute max LDL or MDL for actual values:

tox_list$chem_data <- tox_list$chem_data %>%
  left_join(select(tox_list$chem_info,
                   CAS,
                   MDL = `Maximum method detection level`,
                   LDL = `Maximum laboratory reporting level`),
            by="CAS") %>%
  rowwise() %>%
  mutate(Value = max(MDL, LDL, na.rm = TRUE)) %>%
  select(SiteID, `Sample Date`, CAS, Value) %>%
  distinct()

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

# Basically...need to swap endpoint and site to take advantage of plot_tox_boxplots code

chemicalSummary$site <- chemicalSummary$endPoint



plot_chemical_boxplots_mod <- function(chemicalSummary, 
                                   manual_remove=NULL,
                                   mean_logic = FALSE,
                                   sum_logic = TRUE,
                                   plot_ND = TRUE,
                                   font_size = NA,
                                   title = NA,
                                   pallette = NA,
                                   hit_threshold = NA){
  
  site <- EAR <- sumEAR <- meanEAR <- groupCol <- nonZero <- ".dplyr"
  chnm <- Class <- meanEAR <- x <- y <- ".dplyr"
  
  cbValues <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                "#DCDA4B","#999999","#00FFFF","#CEA226","#CC79A7","#4E26CE",
                "#FFFF00","#78C15A","#79AEAE","#FF0000","#00FF00","#B1611D",
                "#FFA500","#F4426e", "#800000", "#808000")
  
  if(!plot_ND){
    chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
  }
  
  if(length(unique(chemicalSummary$Class)) > length(cbValues)){
    n <- length(unique(chemicalSummary$Class))
    
    if(n > 20 & n<30){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n = 8, name = "Set2"),
                    brewer.pal(n = max(c(3,n-20)), name = "Set1"))
    } else if (n <= 20){
      cbValues <- c(brewer.pal(n = 12, name = "Set3"),
                    brewer.pal(n =  max(c(3,n-12)), name = "Set2"))     
    } else {
      cbValues <- colorRampPalette(brewer.pal(11,"Spectral"))(n)
    }
    
  }
  
  single_site <- length(unique(chemicalSummary$site)) == 1
  
  y_label <- expression(EAR[Chem]~per~ToxCast~Assay)
    
  graphData <- toxEval:::graph_chem_data(chemical_summary = chemicalSummary, 
                               manual_remove=manual_remove,
                               mean_logic=mean_logic,
                               sum_logic=sum_logic)
    
  pretty_logs_new <-  toxEval:::prettyLogs(graphData$meanEAR)
    
  countNonZero <- graphData %>%
    select(chnm, Class, meanEAR) %>%
    group_by(chnm, Class) %>%
    summarize(nonZero = as.character(sum(meanEAR>0)),
              hits = as.character(sum(meanEAR > hit_threshold)))
  
  countNonZero$hits[countNonZero$hits == "0"] <- ""
  
  label <- "# Assays"
  toxPlot_All <- ggplot(data=graphData) 
  
  if(!all(is.na(pallette))){
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=meanEAR, fill=chnm),lwd=0.1,outlier.size=1) +
      scale_fill_manual(values = pallette) +
      theme(legend.position = "none")
  } else {
    toxPlot_All <- toxPlot_All +
      geom_boxplot(aes(x=chnm, y=meanEAR, fill=Class),
                   lwd=0.1,outlier.size=1)
  }

  toxPlot_All <- toxPlot_All +
    scale_y_log10(y_label, labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)  +
    theme_bw() +
    scale_x_discrete(drop = TRUE) +
    geom_hline(yintercept = hit_threshold, linetype="dashed", color="black") +
    theme(axis.text = element_text( color = "black"),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          strip.text.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5))  
  
  if(all(is.na(pallette))){
    toxPlot_All <- toxPlot_All +
      scale_fill_manual(values = cbValues, drop=FALSE) +
      guides(fill=guide_legend(ncol=6)) +
      theme(legend.position="bottom",
            legend.justification = "left",
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
            legend.title=element_blank(),
            legend.key.height = unit(1,"line"))
  }
  
  if(!is.na(font_size)){
    toxPlot_All <- toxPlot_All +
      theme(axis.text = element_text(size = font_size),
            axis.title =   element_text(size=font_size))
  }
  
  #Saving for later!!!!
  if(packageVersion("ggplot2") >= '2.2.1.9000'){
    toxPlot_All <- toxPlot_All +
      coord_flip(clip = "off")
  } else {
    toxPlot_All <- toxPlot_All +
      coord_flip()      
  }
  
  plot_info <- ggplot_build(toxPlot_All)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
    ymax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
  } else {
    ymin <- 10^(layout_stuff$panel_ranges[[1]]$x.range[1])
    ymax <- 10^(layout_stuff$panel_ranges[[1]]$x.range[2])
  }
  
  toxPlot_All_withLabels <- toxPlot_All +
    geom_text(data=countNonZero, aes(x=chnm,label=nonZero, y=ymin), size = ifelse(is.na(font_size),2,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=ymin, label = label, stringsAsFactors = FALSE), 
              aes(x=x,  y=y, label = label),vjust=-0.5,
              size=ifelse(is.na(font_size),3,0.30*font_size)) 
  
  nHitsEP <- countNonZero$hits
  
  if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=countNonZero, aes(x=chnm, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  if(!all(is.na(title))){
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      ggtitle(title)
    
    if(!is.na(font_size)){
      toxPlot_All_withLabels <- toxPlot_All_withLabels +
        theme(plot.title = element_text(size=font_size))
    }
  }
  
  if(!is.na(hit_threshold)) {
    toxPlot_All_withLabels <- toxPlot_All_withLabels +
      geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  return(toxPlot_All_withLabels)
  
}

plot_DL <- plot_chemical_boxplots_mod(chemicalSummary, 
                             font_size = 7, title = " ")

library(gridExtra)

plot_DL_w_cap <- plot_DL +
  labs(caption = bquote(atop(bold("Figure SI-2:") ~ "Exposure activity ratios" ~ (EAR[Chem]) ~ "using ToxCast endpoints and the detection level of chemicals monitored in Great Lakes tributaries, 2010-2013.",
                             "Boxes, 25th to 75th percentiles; dark line, median; whiskers, data within 1.5 x the interquartile range (IQR); circles, values outside 1.5 x the IQR."))) +
  theme(plot.caption = element_text(hjust = 0, size = 6))

dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(plot_DL_w_cap, filename = "plots/SI2_detection_levels1.pdf", width = 9, height = 11)

# line_1 <- expression(bold("Figure SI-2:")*" Exposure activity ratios ("*EAR[Chem]*") using ToxCast endpoints and the detection level of chemicals monitored in Great Lakes")
# line_2 <- expression("tributaries, 2010-2013. Boxes, 25th to 75th percentiles; dark line, median; whiskers, data within 1.5 x the interquartile range (IQR); circles, values outside 1.5 x the IQR.")

# plot_DL_w_cap <- plot_DL +
#   annotation_custom(grid::textGrob(line_1,
#                                    gp = grid::gpar(fontsize = 7)),
#                                    xmin = -6.5, xmax = -6.5, ymin = -4, ymax = -4) +
#   annotation_custom(grid::textGrob(line_2,
#                                    gp = grid::gpar(fontsize = 7)),
#                     xmin = -7.5, xmax = -7.5, ymin = -4, ymax = -4)






# Determine a few things for the text in the manuscript
range(chemicalSummary$EAR)
chemicalSummary[which.max(chemicalSummary$EAR),"chnm"]
chemicalSummary[which.max(chemicalSummary$EAR),"endPoint"]
