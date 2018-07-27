library(toxEval)
library(dplyr)
library(tidyr)
library(readxl)

source(file = "data_setup.R")
AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")


plot_heat_AOPs <- function(chemical_summary,AOP_info,
                                chem_site,
                                mean_logic,
                                sum_logic){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  graphData <- toxEval:::graph_chem_data(chemical_summary, 
                               mean_logic=mean_logic,
                               sum_logic = sum_logic)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  graphData <- chemical_summary %>%
    select(-Class) %>%
    left_join(select(AOP_info, AOP, endPoint=`Endpoint(s)`, Class=`X__1`), by="endPoint") %>%
    group_by(site,date,AOP, Class) %>%
    summarise(sumEAR=sum(EAR)) %>%
    data.frame() %>%
    group_by(site, AOP, Class) %>%
    summarise(meanEAR=ifelse(mean_logic,mean(sumEAR),max(sumEAR))) %>%
    data.frame() 
  
  graphData <- graphData %>%
    left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
              by=c("site"="SiteID"))
  
  graphData$Class[is.na(graphData$Class)] <- "AOP not defined"
  graphData$AOP[is.na(graphData$AOP)] <- "None"
  
  class_order <- graphData %>%
    filter(!(Class %in% c("AOP not defined","Not environmentally relevant"))) %>%
    group_by(Class) %>%
    summarize(maxEAR = max(meanEAR)) %>%
    arrange(desc(maxEAR)) %>%
    pull(Class)
  
  graphData$AOP <- factor(graphData$AOP)
  graphData$Class <- factor(graphData$Class, levels = c(class_order,"AOP not defined","Not environmentally relevant"))
  
  fill_text <- ifelse(mean_logic, "mean", "max")
  fill_text <- bquote(italic(.(fill_text))~group("[", EAR["[" * j * "]"] , "]"))
  
  heat <- ggplot(data = graphData) +
    geom_tile(aes(x = `Short Name`, y=AOP, fill=meanEAR, color="")) +
    theme_bw() +
    theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) +
    ylab("AOP ID") +
    xlab("") +
    labs(fill=fill_text) +
    scale_color_manual(values=NA) +
    scale_fill_gradient( trans = 'log',
                         low = "white", high = "steelblue",
                         breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,5),
                         na.value = 'khaki',labels=toxEval:::fancyNumbers2) +
    facet_grid(Class ~ site_grouping, scales="free", space="free") +
    labs(caption = bquote(atop(bold("Figure SI-6:") ~ "Maximum exposure activity ratios (" * .(fill_text) *") for each adverse outcome pathway (AOP) identified", 
                               "from evaluation of chemistry data at monitored Great Lakes tributaries, 2010-2013 ("*italic("j = samples") *").                          "))) +
    theme(strip.text.y = element_text(angle=0, hjust=0), 
          strip.background = element_rect(fill="transparent", colour = NA),
          # axis.text.y = element_text(face=ifelse(levels(graphData$category) %in% c("Total"),"bold","italic")),
          panel.spacing = unit(0.05, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA),
          plot.caption = element_text(hjust = 0)) +
    guides(color=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
           fill = guide_colorbar(order=1))
  
  return(heat)
  
}

aop_heat <- plot_heat_AOPs(chemicalSummary, AOP_info, tox_list$chem_site, 
               sum_logic = FALSE, mean_logic = FALSE)

ggsave(aop_heat, file="plots/SI6_AOP_heat.pdf", height = 9, width = 11)
