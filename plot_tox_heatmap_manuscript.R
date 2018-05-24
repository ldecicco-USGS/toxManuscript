
plot_tox_heatmap_manuscript <- function(chemicalSummary, 
                             chem_site, 
                             category = "Biological",
                             breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                             manual_remove = NULL,
                             mean_logic = FALSE,
                             sum_logic = TRUE,
                             plot_ND = TRUE, 
                             font_size = NA,
                             title = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))
  
  

  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  fill_label <- ifelse(mean_logic, "Mean EAR", "Max EAR")
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  if(!plot_ND){
    chemicalSummary <- chemicalSummary[chemicalSummary$EAR > 0,]
  }
  
  if(category == "Chemical"){
    plot_back <- plot_heat_chemicals_manuscript(chemicalSummary=chemicalSummary, 
                                     mean_logic=mean_logic,
                                     sum_logic=sum_logic,
                                     chem_site=chem_site)
    
  } else {
    
    graphData <- tox_boxplot_data(chemicalSummary = chemicalSummary,
                                  category = category,
                                  manual_remove = manual_remove,
                                  mean_logic = mean_logic,
                                  sum_logic = sum_logic)
    
    graphData <- graphData %>%
      left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
                by=c("site"="SiteID"))
    
    
    plot_back <- ggplot(data = graphData) +
      geom_tile(aes(x = `Short Name`, y=category, fill=meanEAR)) +
      theme_bw() +
      theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 0.975)) +
      ylab("") +
      xlab("") +
      labs(fill=fill_label) +
      scale_fill_gradient( guide = "legend",
                           trans = 'log',
                           low = "white", high = "steelblue",
                           breaks=breaks,
                           na.value = 'transparent',labels=fancyNumbers2) +
      facet_grid(. ~ site_grouping, scales="free", space="free") +
      theme(strip.text.y = element_text(angle=0, hjust=0), 
            strip.background = element_rect(fill="transparent", colour = NA),
            panel.spacing = unit(0.05, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA))
    
  }
  
  if(!is.na(font_size)){
    plot_back <- plot_back +
      theme(axis.text = element_text(size = font_size),
            strip.text = element_text(size = font_size))
  }
  
  if(!is.na(title)){
    plot_back <- plot_back +
      ggtitle(title)
    
    if(!is.na(font_size)){
      plot_back <- plot_back +
        theme(plot.title = element_text(size=font_size))
    }
  }  
  
  return(plot_back)
}


plot_heat_chemicals_manuscript <- function(chemicalSummary,
                                chem_site,
                                mean_logic,
                                sum_logic){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  priority_chems <- chemicalSummary %>%
    group_by(site,date,chnm) %>%
    summarise(sumEAR = sum(EAR)) %>%
    group_by(chnm,site) %>%
    summarise(maxEAR = max(sumEAR, na.rm = TRUE)) %>%
    filter(maxEAR > ear_thresh) %>%
    group_by(chnm) %>%
    summarise(sites = n_distinct(site)) %>%
    filter(sites > siteThres) %>%
    pull(chnm)
  

  graphData <- graph_chem_data(chemicalSummary, 
                               mean_logic=mean_logic,
                               sum_logic = sum_logic)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  graphData <- graphData %>%
    left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
              by=c("site"="SiteID"))
  
  fill_text <- ifelse(mean_logic, "Mean EAR", "Max EAR")
  
  # colors <- ifelse(levels(graphData$chnm) %in% as.character(priority_chems),"red","black")
  graphData$v <- ifelse(graphData$chnm %in% as.character(priority_chems), "red", "black" )
  
  heat <- ggplot(data = graphData) +
    geom_tile(aes(x = `Short Name`, y=chnm, fill=meanEAR)) +
    theme_bw() +
    ylab("") +
    xlab("") +
    labs(fill=fill_text) +
    scale_fill_gradient( guide = "legend",
                         trans = 'log',
                         low = "white", high = "steelblue",
                         breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                         na.value = 'transparent',labels=toxEval:::fancyNumbers2) +
    facet_grid(Class ~ site_grouping, scales="free", space="free") +
    theme(strip.text.y = element_text(angle=0, hjust=0), 
          strip.background = element_rect(fill="transparent", colour = NA),
          axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1),
          # axis.text.y = element_text(color = c(rep("black",2),"red",rep("black",2),"blue")),
          panel.spacing = unit(0.05, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA)) 
    # expand_limits(x = -10) +
    # geom_rect(ymin = 0, ymax = length(levels(graphData$chnm))+1, xmin = -5, xmax = 0, fill = "white") +
    # # # fake axis layer, aligned below y = 0
    # geom_text(data = distinct(filter(select(graphData, v, chnm, site_grouping, Class), site_grouping == "Lake Superior")), aes(color = v, label = chnm, y=chnm), x = 0) +
    # # # # specify the font colours for fake axis
    # scale_colour_manual(values = c("black", "red"), guide = F) +
    # # # hide the actual x-axis text / ticks
    # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    # # 
    # 
  return(heat)
  
}



