source(file = "data_setup.R")

######################################################3
AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

chemicalSummary <- left_join(chemicalSummary, AOP, by="endPoint") %>%
  select(-Bio_category) %>%
  rename(Bio_category = ID)

plot_aop <- function(chemicalSummary, 
                     category = "AOP",
                     filterBy = "All",
                     manual_remove = NULL,
                     hit_threshold = 0.1,
                     mean_logic = FALSE, 
                     font_size = NA,
                     title = NA,
                     pallette = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical","AOP"))
  
  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else if (category == "AOP"){
    chemicalSummary$category <- as.character(chemicalSummary$ID)
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
  
  graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                                category = category,
                                manual_remove = manual_remove,
                                mean_logic = mean_logic)
  
  countNonZero <- graphData %>%
    group_by(category) %>%
    summarise(nonZero = as.character(length(unique(site[meanEAR>0])))) %>%
    data.frame() 
  
  label <- "# Sites"
  
  bioPlot <- ggplot(data = graphData)+
    coord_flip() +
    theme_bw() +
    xlab("AOP ID") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          axis.text.y = element_text(color = "black", vjust = 0.2), 
          axis.text.x = element_text(color = "black", vjust = 0, margin = margin(-0.5,0,0,0)),
          panel.border = element_blank(),
          axis.ticks = element_blank()) +  
    scale_y_log10("Maximum EAR Per Site",labels=toxEval:::fancyNumbers)
  
  if(!all(is.na(pallette))){
    bioPlot <- bioPlot +
      geom_boxplot(aes(x=category, y=meanEAR, fill = category),lwd=0.1,outlier.size=1) +
      scale_fill_manual(values = cbValues) +
      theme(legend.position = "none")
  } else {
    bioPlot <- bioPlot +
      geom_boxplot(aes(x=category, y=meanEAR),lwd=0.1,outlier.size=1, fill = "steelblue")      
  }
  
  if(!is.na(font_size)){
    bioPlot <- bioPlot +
      theme(axis.text = element_text(size = font_size))
  }
  
  plot_info <- ggplot_build(bioPlot)
  layout_stuff <- plot_info$layout
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    xmin <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[1])
    xmax <- 10^(layout_stuff$panel_scales_y[[1]]$range$range[2])
    ymax <- layout_stuff$panel_scales_x[[1]]$range$range[2]
  } else {
    xmin <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[1]))
    xmax <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[2]))
    ymax <- suppressWarnings(layout_stuff$panel_ranges[[1]]$y.range[2])
  }
  
  bioPlot_w_labels <- bioPlot + 
    geom_text(data=countNonZero, aes(x=category, y=xmin,label=nonZero),size=ifelse(is.na(font_size),3,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=xmin, label = label, stringsAsFactors = FALSE), 
              aes(x = x,  y=y, label = label),
              size=ifelse(is.na(font_size),3,0.30*font_size)) 
  
  if(!is.na(title)){
    bioPlot_w_labels <- bioPlot_w_labels +
      ggtitle(title)
    
    if(!is.na(font_size)){
      bioPlot_w_labels <- bioPlot_w_labels +
        theme(plot.title = element_text(size=font_size))
    }
  }  
  
  return(bioPlot_w_labels)
  
}

chem_1 <- "Atrazine"

chem_sum_1 <- filter(chemicalSummary, chnm == chem_1) %>%
  select(-Bio_category) %>%
  rename(Bio_category = ID)

aop_graph <- plot_aop(chem_sum_1, chem_1, category = "Biological", title = chem_1)

gb <- ggplot2::ggplot_build(aop_graph)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel"] <- "off"

png("atrazine_AOP.png", width = 600, height = 600, res = 142)
grid::grid.draw(gt)
dev.off()

