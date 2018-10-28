library(toxEval)
library(dplyr)
library(tidyr)

source(file = "data_setup.R")

EARthresh <- 0.001
siteThresh <- 10

SiteChemsCount <- filter(chemicalSummary,EAR>0) %>%
  group_by(shortName,date,chnm) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(shortName,chnm) %>%
  summarize(EARsummax = max(EARsum)) %>%
  filter(EARsummax > EARthresh) %>%
  group_by(chnm) %>%
  summarize(numSites = n_distinct(shortName)) %>%
  filter(numSites >= siteThresh)

priorityChems <- as.character(SiteChemsCount$chnm)
priorityChemRows <- which(chemicalSummary$chnm %in% priorityChems)

plot_heat_chemicals_man <- function(chemical_summary,
                                chem_site,priority_chems,
                                mean_logic,
                                sum_logic,
                                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10)){
  
  SiteID <- site_grouping <- `Short Name` <- chnm <- maxEAR <- ".dplyr"
  site <- EAR <- sumEAR <- meanEAR <- ".dplyr"
  
  graphData <- toxEval:::graph_chem_data(chemical_summary, 
                               mean_logic=mean_logic,
                               sum_logic = sum_logic)
  
  if(!("site_grouping" %in% names(chem_site))){
    chem_site$site_grouping <- "Sites"
  }
  single_site <- length(unique(chemical_summary$site)) == 1
  
  fill_text <- ifelse(mean_logic, "mean", "max")
  fill_text <- bquote(italic(.(fill_text))* EAR[ChemSite])
  
  graphData <- graphData %>%
    left_join(chem_site[, c("SiteID", "site_grouping", "Short Name")],
              by=c("site"="SiteID"))
  
  # This requires non-detects to be 0. If that changes we'll need to update:
  graphData$meanEAR[graphData$meanEAR == 0] <- NA
  
  complete_data_filled <- get_complete_set(chemical_summary, graphData, chem_site)
  
  any_missing <- nrow(complete_data_filled) > nrow(graphData)
  any_non_detects <- any(is.na(graphData$meanEAR))
  
  color_df <- data.frame(chnm = levels(graphData$chnm)) %>%
    mutate(color = ifelse(as.character(chnm) %in% as.character(priority_chems), "red", "black")) %>%
    left_join(distinct(select(graphData, chnm, Class)), by="chnm") %>%
    mutate(site_grouping = factor(rep("Lake Superior",length(levels(graphData$chnm))), levels = levels(graphData$site_grouping)),
           chnm = factor(chnm, levels = levels(graphData$chnm)))
  
  heat <- ggplot(data = graphData) +
    geom_point(data = complete_data_filled, aes(x = `Short Name`, y=chnm, shape=""), size = 1 ) +
    geom_tile(aes(x = `Short Name`, y=chnm, fill=meanEAR, alpha = "")) +
    theme_bw() +
    theme(axis.text.x = element_text( angle = 90,vjust=0.5,hjust = 1)) +
    ylab("") +
    xlab("") +
    labs(fill = fill_text, caption = bquote(italic("j = samples"))) +
    scale_fill_gradient( na.value = 'khaki',
                         trans = 'log', low = "white", high = "steelblue",
                         breaks=breaks,
                         labels=toxEval:::fancyNumbers2, guide = "colourbar") +
    scale_alpha_manual(values=NA) +
    scale_shape_manual(values=4) +
    facet_grid(Class ~ site_grouping, scales="free", space="free") +
    theme(strip.text.y = element_text(angle=0, hjust=0), 
          strip.background = element_rect(fill="transparent", colour = NA),
          # axis.text.y = element_text(face=ifelse(levels(graphData$category) %in% c("Total"),"bold","italic")),
          panel.spacing = unit(0.05, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = "transparent")) +
    coord_cartesian(clip = "off") +
    # # # fake axis layer, aligned below y = 0
    geom_text(data = color_df, 
              aes(color = color, label = chnm, y=chnm), x = 0.2,
              size = 2.5, # match font size to theme
              hjust = 1, vjust = 0.3) +
    # # # # specify the font colours for fake axis
    scale_colour_manual(values = c("black", "red"), guide = F) +
    # # # hide the actual x-axis text / ticks
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_text(margin = margin(r = 125)))
  
  if(any_non_detects & any_missing){
    heat <- heat +
      guides(colour=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
             shape=guide_legend("Missing", order = 3),
             fill = guide_colorbar(order=1)) 
  } else if (any_non_detects){
    heat <- heat +
      guides(alpha=guide_legend("Non-detects", override.aes=list(colour="khaki", fill="khaki"), order = 2),
             fill = guide_colorbar(order=1),
             shape = "none")    
  } else if (any_missing){
    heat <- heat +
      guides(shape=guide_legend("Missing", order = 2),
             fill = guide_colorbar(order=1),
             alpha = "none")
  } else {
    heat <- heat +
      guides(fill = guide_colorbar(order=1),
             shape = "none",
             alpha = "none")
  }
  
  return(heat)
  
}

# There's probably a faster way to do this:
get_complete_set <- function(chemical_summary, graphData, chem_site){
  
  complete_data <- select(chem_site, `Short Name`, site_grouping)
  complete_data_filled <- data.frame()
  
  for(chms in levels(chemical_summary$chnm)){
    complete_data$chnm <- chms
    complete_data_filled <- dplyr::bind_rows(complete_data_filled, complete_data)
  }
  complete_data_filled$chnm <- factor(complete_data_filled$chnm, levels = levels(graphData$chnm))
  complete_data_filled <- left_join(complete_data_filled, distinct(select(chemical_summary, chnm, Class)), by="chnm")
  return(complete_data_filled)
}

heat_map <- plot_heat_chemicals_man(chemicalSummary, tox_list$chem_site, priorityChems,
                                    mean_logic = FALSE, sum_logic = TRUE)

heat_map_w_cap <- heat_map +
  labs(caption = bquote(atop(bold("Figure SI-4:") ~ "Maximum exposure acivity ratios (" * 
                               italic("max")~(EAR[ChemSite]) *
                          ") at monitored Great Lakes tributaries, 2010-2013. Chemicals in red ",
                          "indicate those that were present at a minimum of 10 sites and had" ~
                          italic("max")~(EAR[ChemSite]) ~
                          ">"~ 10^-3 ~" at a minimum of 5 sites.        "))) +
  theme(plot.caption = element_text(hjust = 0.85))

ggsave(heat_map_w_cap, filename = "plots/SI4_heat_map.pdf", height = 9, width = 11)

