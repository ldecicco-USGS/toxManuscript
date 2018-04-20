source(file = "data_setup.R")

######################################################3
AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

chemicalSummary <- left_join(chemicalSummary, AOP, by = "endPoint")

chemicalSummary$has_AOP <- "Associated AOP"
chemicalSummary$has_AOP[is.na(chemicalSummary$ID)] <- "Undefined AOP"


plot_tox_endpoints_special <- function(chemicalSummary, 
                                       category = "Biological",
                                       filterBy = "All",
                                       manual_remove = NULL,
                                       hit_threshold = 0.1,
                                       mean_logic = FALSE, 
                                       font_size = NA,
                                       title = NA,
                                       pallette = NA,
                                       include_legend = TRUE,
                                       include_axis = TRUE){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
  
  if(filterBy != "All"){
    if(!(filterBy %in% unique(chemicalSummary$category))){
      stop("filterBy argument doesn't match data")
    }
    
    chemicalSummary <- chemicalSummary %>%
      filter_(paste0("category == '", filterBy,"'"))
  }
  
  graphData <- chemicalSummary %>%
    group_by(site,date,category,endPoint,has_AOP) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, category,endPoint,has_AOP) %>%
    summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
    data.frame() %>%
    mutate(category=as.character(category))
  
  countNonZero <- graphData %>%
    group_by(endPoint) %>%
    summarise(nonZero = as.character(sum(meanEAR>0)),
              hits = as.character(sum(meanEAR > hit_threshold)))
  
  countNonZero$hits[countNonZero$hits == "0"] <- ""
  
  namesToPlotEP <- as.character(countNonZero$endPoint)
  nSamplesEP <- countNonZero$nonZero
  nHitsEP <- countNonZero$hits
  
  orderColsBy <- graphData %>%
    group_by(endPoint) %>%
    summarise(median = quantile(meanEAR[meanEAR != 0],0.5)) %>%
    arrange(median)
  
  orderedLevelsEP <- orderColsBy$endPoint
  
  if(any(is.na(orderColsBy$median))){
    orderedLevelsEP <- c(orderColsBy$endPoint[is.na(orderColsBy$median)],
                         orderColsBy$endPoint[!is.na(orderColsBy$median)])
  }
  
  graphData$endPoint <- factor(graphData$endPoint, levels = orderedLevelsEP)
  
  stackedPlot <- ggplot(graphData)+
    scale_y_log10("Maximum EAR Per Site",labels=toxEval:::fancyNumbers) +
    theme_minimal() +
    xlab("") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  
  if(!include_axis){
    stackedPlot <- stackedPlot +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  } else {
    stackedPlot <- stackedPlot +
      theme(axis.text.y = element_text(color = "black", vjust = 0.2), 
            axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)))
  }
  
  if(!is.na(hit_threshold)){
    stackedPlot <- stackedPlot +
      geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
  }
  
  
  if(!all(is.na(pallette))){
    stackedPlot <- stackedPlot +
      geom_boxplot(aes(x=endPoint, y=meanEAR, fill = has_AOP, color=has_AOP),outlier.size = 0.1 ) +
      scale_fill_manual(values = pallette) +
      scale_color_manual(values = pallette) 
    
    if(include_legend){
      stackedPlot <- stackedPlot +
        theme(legend.title=element_blank())
    } else {
      stackedPlot <- stackedPlot +
        theme(legend.position = "none")
    }
    
  } else {
    stackedPlot <- stackedPlot +
      geom_boxplot(aes(x=endPoint, y=meanEAR), fill = "steelblue") 
  }
  
  if(!is.na(font_size)){
    stackedPlot <- stackedPlot +
      theme(axis.text = element_text(size = font_size))
  }
  
  if(!is.na(title)){
    stackedPlot <- stackedPlot +
      ggtitle(title)
    
    if(!is.na(font_size)){
      stackedPlot <- stackedPlot +
        theme(plot.title = element_text(size=font_size))
    }
  } 
  
  stackedPlot <- stackedPlot +
    coord_flip()
  
  return(stackedPlot)
}

pdf("AOP_endpoint_better.pdf", width = 9, height = 11)
for(chemical in rev(levels(chemicalSummary$chnm))){
  
  chem_sum_1 <- filter(chemicalSummary, chnm == chemical)
  ep_graph2 <- plot_tox_endpoints_special(chem_sum_1, 
                                          category = "Chemical",
                                          title=chemical,
                                          hit_threshold = NA,
                                          include_axis = TRUE,
                                          include_legend = TRUE,
                                          pallette = c("#e69f00","#56B4E9"))
  print(ep_graph2)
  
}
dev.off()

chem_1 <- "Atrazine"

chem_sum_1 <- filter(chemicalSummary, chnm == chem_1)

at_plot <- plot_tox_endpoints_special(chem_sum_1, 
                                      category = "Chemical",
                                      title=chem_1,
                                      hit_threshold = NA,
                                      include_axis = TRUE,
                                      include_legend = TRUE,
                                      pallette = c("#e69f00","#56B4E9"))

ggsave(at_plot, filename = "atrazine.png", width = 10, height = 5)

at_plot <- plot_tox_endpoints_special(chem_sum_1, 
                                      category = "Chemical",
                                      title=chem_1,
                                      hit_threshold = NA,
                                      include_axis = TRUE,
                                      include_legend = FALSE,
                                      pallette = c("#e69f00","#56B4E9"))
ggsave(at_plot, filename = "atrazine_noLegend.png", width = 8, height = 5)

chem_2 <- "Metolachlor"
chem_sum_2 <- filter(chemicalSummary, chnm == chem_2)

met_plot <- plot_tox_endpoints_special(chem_sum_2, 
                                       category = "Chemical",
                                       title=chem_2,
                                       hit_threshold = NA,
                                       include_axis = TRUE,
                                       include_legend = TRUE,
                                       pallette = c("#e69f00","#56B4E9"))

ggsave(met_plot, filename = "metolachlor.png", width = 10, height = 5)

met_plot <- plot_tox_endpoints_special(chem_sum_2, 
                                       category = "Chemical",
                                       title=chem_2,
                                       hit_threshold = NA,
                                       include_axis = TRUE,
                                       include_legend = FALSE,
                                       pallette = c("#e69f00","#56B4E9"))
ggsave(met_plot, filename = "metolachlor_noLegend.png", width = 8, height = 5)


#####################################
# Big graph:

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  filter(endPoint %in% unique(chemicalSummary$endPoint)) %>%
  select(endPoint) %>%
  distinct()

endPoints_with_AOP <- AOP$endPoint
endPoints_without <- chemicalSummary$endPoint[!(chemicalSummary$endPoint %in% endPoints_with_AOP)]

chemicalSummary$has_AOP <- "Associated AOP"
chemicalSummary$has_AOP[chemicalSummary$endPoint %in% endPoints_without] <- "Undefined AOP"

big_ep <- plot_tox_endpoints_special(chemicalSummary, 
                                     category = "Chemical",
                                     include_axis  = FALSE,
                                     include_legend = FALSE,
                                     hit_threshold = NA, 
                                     pallette = c("#e69f00","#56B4E9"))

ggsave(big_ep, filename = "Big_endpoint_no_legend.png")



