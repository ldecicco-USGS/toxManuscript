#' EndPoint boxplots
#' 
#' This function creates a set of boxplots for each endPoint based on a selected 
#' filter. The data is first filtered down to what group is specified in the 'filterBy' 
#' argument. The 'filterBy' argument must match one of the unique options in the 'category'. 
#' So for example, if the 'category' is "Chemical Class", then the 'filterBy' argument must 
#' be one of the defined "Chemical Class" options such as "Herbicide". 
#'
#' After the data is filtered, a boxplot is generated for each endPoint. The EAR values 
#' that are used to create the boxplots are still the mean or max (as defined by 'mean_logic') 
#' of each site.
#' 
#' @param chemicalSummary data frame from \code{get_chemical_summary}
#' @param category either "Biological", "Chemical Class", or "Chemical"
#' @param filterBy character either "All" or one of the filtered categories.
#' @param manual_remove vector of categories to remove
#' @param mean_logic character. Options are "mean", "max", or "noSum". 
#' TRUE will default to "mean" and FALSE to "max". The default value is "mean". 
#' The most appropriate use of "noSum" is for non-ToxCast benchmarks. In this case
#' the values plotted are the overall max of the sample (not the max of the sum
#' of the sample).
#' @param hit_threshold numeric threshold defining a "hit"
#' @param font_size numeric to adjust the axis font size
#' @param title character title for plot. 
#' @param pallette vector of color pallette for fill. Can be a named vector
#' to specify specific color for specific categories. 
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @importFrom utils packageVersion
#' @importFrom dplyr full_join filter mutate select left_join right_join
#' @examples
#' # This is the example workflow:
#' path_to_tox <-  system.file("extdata", package="toxEval")
#' file_name <- "OWC_data_fromSup.xlsx"
#'
#' full_path <- file.path(path_to_tox, file_name)
#' 
#' \dontrun{
#' tox_list <- create_toxEval(full_path)
#' ACClong <- get_ACC(tox_list$chem_info$CAS)
#' ACClong <- remove_flags(ACClong)
#' 
#' cleaned_ep <- clean_endPoint_info(endPointInfo)
#' filtered_ep <- filter_groups(cleaned_ep)
#' chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
#'
#' plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")
#' plot_tox_endpoints(chemicalSummary, category = "Chemical Class", filterBy = "PAHs")
#' plot_tox_endpoints(chemicalSummary, category = "Chemical", filterBy = "Atrazine")
#' 
#'  # To turn off clipping:
#' class_plot <- plot_tox_endpoints(chemicalSummary, filterBy = "Cell Cycle")
#' gb <- ggplot2::ggplot_build(class_plot)
#' gt <- ggplot2::ggplot_gtable(gb)
#' 
#' gt$layout$clip[gt$layout$name=="panel"] <- "off"
#' 
#' grid::grid.draw(gt) 
#' }
plot_tox_endpoints_manuscript <- function(chemicalSummary, 
                              category = "Biological",
                              filterBy = "All",
                              manual_remove = NULL,
                              hit_threshold = NA,
                              mean_logic = FALSE, 
                              sum_logic = TRUE,
                              font_size = NA,
                              title = NA,
                              pallette = NA){
  
  match.arg(category, c("Biological","Chemical Class","Chemical"))

  site <- endPoint <- EAR <- sumEAR <- meanEAR <- x <- y <- ".dplyr"
  
  if(category == "Biological"){
    chemicalSummary$category <- chemicalSummary$Bio_category
  } else if(category == "Chemical Class") {
    chemicalSummary$category <- chemicalSummary$Class
  } else {
    chemicalSummary$category <- chemicalSummary$chnm
  }
      
  single_site <- length(unique(chemicalSummary$site)) == 1
  
  y_label <- toxEval:::fancyLabels(category, mean_logic, sum_logic, single_site)
  
  if(filterBy != "All"){
    if(!(filterBy %in% unique(chemicalSummary$category))){
      stop("filterBy argument doesn't match data")
    }
    
    chemicalSummary <- chemicalSummary %>%
      filter_(paste0("category == '", filterBy,"'"))
  }

  graphData <- chemicalSummary %>%
    group_by(site,date,endPoint,has_AOP) %>%
    summarise(sumEAR=sum(EAR)) %>%
    group_by(site, endPoint,has_AOP) %>%
    summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) 
    
  
  pretty_logs_new <- toxEval:::prettyLogs(graphData$meanEAR)
    
  chem_ns <- chemicalSummary %>%
    group_by(endPoint) %>%
    summarize(nChems = length(unique(chnm)))
  
  nChemsEP <- chem_ns$nChems  

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
    scale_y_log10(expression(EAR[mixture]),labels=toxEval:::fancyNumbers,breaks=pretty_logs_new) +
    theme_minimal() +
    xlab("ToxCast Assay Name") +
    theme(axis.text.y = element_text(vjust = .25,hjust=1)) +
    theme(legend.title = element_blank()) +
    geom_boxplot(aes(x=endPoint, y=meanEAR, fill = has_AOP)) 
  
  if(!is.na(hit_threshold)){
    stackedPlot <- stackedPlot +
      geom_hline(yintercept = hit_threshold, linetype="dashed", color="black")
  }
  
  if(!all(is.na(pallette))){
    stackedPlot <- stackedPlot +
      scale_fill_manual(values = pallette) 
  } 
    
  plot_layout <- ggplot_build(stackedPlot)$layout    
  
  if(packageVersion("ggplot2") >= "2.2.1.9000"){
    ymin <- 10^(plot_layout$panel_scales_y[[1]]$range$range[1])
    ymax <- 10^(plot_layout$panel_scales_y[[1]]$range$range[2])
    
    xmax <- plot_layout$plot_layout$panel_scales_x[[1]]$range$range[1]
    xmin <- plot_layout$plot_layout$panel_scales_x[[1]]$range$range[2]
  } else {
    ymin <- 10^(plot_layout$panel_ranges[[1]]$y.range)[1]
    ymax <- 10^(plot_layout$panel_ranges[[1]]$y.range)[2]

    xmax <- plot_layout$panel_ranges[[1]]$x.range[2]
    xmin <- plot_layout$panel_ranges[[1]]$x.range[1]
  }
  
  label <- c("# Sites","# Chemicals")

  stackedPlot <- stackedPlot +
    geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymin,label=nSamplesEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
    geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nChemsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
    geom_text(data=data.frame(x = Inf, y=c(ymin,ymax), label = label, stringsAsFactors = FALSE), 
              aes(x = x,  y=y, label = label), vjust = -0.5,
              size=ifelse(is.na(font_size),3,0.30*font_size))     
  
  if(isTRUE(sum(as.numeric(nHitsEP), na.rm = TRUE) > 0)) {
    stackedPlot <- stackedPlot +
      geom_text(data=data.frame(), aes(x=namesToPlotEP, y=ymax,label=nHitsEP),size=ifelse(is.na(font_size),3,0.30*font_size)) +
      geom_text(data=data.frame(x = Inf, y=ymax, label = "# Hits", stringsAsFactors = FALSE), 
              aes(x = x,  y=y, label = label),
              size=ifelse(is.na(font_size),3,0.30*font_size))
  }
  
  if(!is.na(hit_threshold)) {
    stackedPlot <- stackedPlot +
      geom_text(data=data.frame(x = Inf, y=hit_threshold, label = "Threshold", stringsAsFactors = FALSE), 
                aes(x = x,  y=y, label = label),
                size=ifelse(is.na(font_size),3,0.30*font_size))
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