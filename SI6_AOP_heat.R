library(toxEval)
library(dplyr)
library(tidyr)
library(readxl)

source(file = "data_setup.R")
AOP_info <- read_xlsx("SI_6_AOP_relevance With Short AOP name.xlsx", sheet = "SI_AOP_relevance")

ear_thresh <- 0.001
siteThres <- 10
# ep_percent_thres <- 0.5

AOP <- read.csv("AOP_crosswalk_Dec_2018.csv", stringsAsFactors = FALSE) %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- data.table::fread("AOP_relevance.csv", data.table = FALSE) %>%
  distinct()

eps_with_ids <- unique(AOP$endPoint)

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= ear_thresh) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThres) %>%
  mutate(hasAOP = endPoint %in% eps_with_ids)

priority_endpoints <- endpoints_sites_hits$endPoint[endpoints_sites_hits$hasAOP]

AOP_full_info <- relevance %>%
  select(-`Endpoint(s)`) %>%
  left_join(AOP, by=c("AOP"="ID")) %>%
  left_join(select(AOP_info, `Abbreviated AOP description` = `...5`, AOP), by="AOP") %>%
  select(AOP, Relevant, Rationale, `Abbreviated AOP description`, `Tox Cast Endpoints`=endPoint) %>%
  arrange(AOP) %>%
  distinct()

AOP_full_info_prior <- AOP_full_info %>%
  filter(`Tox Cast Endpoints` %in% priority_endpoints )

plot_heat_AOPs <- function(chemical_summary,AOP_info,
                                chem_site,
                                mean_logic,
                                sum_logic){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"

  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  
  graphData <- chemical_summary %>%
    select(-Class) %>%
    left_join(select(AOP_info, AOP, 
                     endPoint=`Tox Cast Endpoints`, 
                     Class=`Abbreviated AOP description`), by="endPoint") %>%
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
  fill_text <- bquote(EAR[SiteAOP])
  
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
    labs(caption = bquote(atop(bold("Figure SI-6:") ~ "Exposure activity ratios (" * .(fill_text) *") for each adverse outcome pathway (AOP) identified", 
                               "from evaluation of chemistry data at monitored Great Lakes tributaries, 2010-2013.                          "))) +
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

aop_heat <- plot_heat_AOPs(chemicalSummary, 
                           AOP_full_info, 
                           tox_list$chem_site, 
               sum_logic = FALSE, mean_logic = FALSE)

# ggsave(aop_heat, file="plots/SI6_AOP_heat.pdf", height = 9, width = 11)

ggsave(aop_heat, file="plots/SI6_AOP_heat_all.pdf", height = 20, width = 11)

aop_heat2 <- plot_heat_AOPs(filter(chemicalSummary, endPoint %in% priority_endpoints), 
                           AOP_full_info_priority, 
                           tox_list$chem_site, 
                           sum_logic = FALSE, mean_logic = FALSE)

ggsave(aop_heat2, file="plots/SI6_AOP_heat_priorityEPs.pdf", height = 9, width = 11)
