library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(dataRetrieval)
library(cowplot)
library(grid)

####################################
source(file = "data_setup.R")

ear_thresh <- 0.001
siteThres <- 10
ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)

relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.) %>%
  mutate(ID = factor(ID))

boxData_max <- chemicalSummary %>%
  left_join(AOP, by="endPoint") %>%
  group_by(ID, chnm, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)]) %>%
  ungroup() %>%
  filter(maxEAR > 0) %>%
  mutate(ID = as.factor(ID))

boxData_tots <- boxData_max %>%
  group_by(ID,site,date) %>%
  summarize(total = sum(maxEAR))  %>%
  ungroup() 

boxData <- boxData_tots %>%
  mutate(ID = as.factor(ID)) %>%
  group_by(ID, site) %>%
  summarize(maxMaxEAR = max(total, na.rm = TRUE),
            date_picked = date[which.max(total)]) %>%
  ungroup() %>%
  filter(!is.na(ID),
         maxMaxEAR > ear_thresh) 

priority_AOPs <- boxData %>%
  group_by(ID) %>%
  summarise(siteDet = n_distinct(site)) %>%
  filter(siteDet >= siteThres)

boxData <- filter(boxData, ID %in% priority_AOPs$ID)
boxData_tots <- filter(boxData_tots, ID %in% priority_AOPs$ID)
boxData_max <- filter(boxData_max, ID %in% priority_AOPs$ID)

relevance$ID <- factor(as.character(relevance$ID), levels = levels(boxData$ID))

boxData <- boxData %>%
  left_join(select(relevance, ID, Relevant), by="ID")

fractions <- boxData_tots %>%
  left_join(boxData_max, by=c("ID","site","date")) %>%
  right_join(boxData, by=c("ID","site","date"="date_picked")) %>%
  group_by(ID, site, date, chnm, endPoint_used) %>%
  summarize(fraction = maxEAR/total) %>%
  ungroup()

endpoints_that_contribute <- fractions %>%
  filter(fraction > ep_percent_thres) %>%
  select(endPoint_used) %>%
  distinct() %>%
  pull(endPoint_used)

chem_sum_AOP <- boxData_max %>%
  right_join(select(boxData, -ID), by=c("site","date"="date_picked")) %>%
  filter(endPoint_used %in% endpoints_that_contribute) %>%
  group_by(ID, endPoint_used) %>%
  summarize(maxEAR = max(maxEAR, na.rm = TRUE),
            meanEAR = mean(maxEAR, na.rm = TRUE),
            medianEAR = median(maxEAR, na.rm = TRUE)) %>%
  data.frame() %>%
  filter(!is.na(ID)) %>%
  mutate(ID = as.factor(ID),
         endPoint = as.factor(endPoint_used)) 

nSites <- boxData %>%
  group_by(ID) %>%
  summarize(sitehits = sum(maxMaxEAR > ear_thresh)) %>%
  filter(!is.na(ID),
         ID %in% priority_AOPs$ID) %>%
  mutate(ID = as.factor(ID))

chem_sum_AOP <- filter(chem_sum_AOP, ID %in% priority_AOPs$ID)
chem_sum_AOP$endPoint <- droplevels(chem_sum_AOP$endPoint )

pretty_logs_new <- toxEval:::prettyLogs(boxData$maxMaxEAR)

y_label <- bquote("max" ~ 
                    group("[", 
                          group("(",
                                sum(" "  ~ EAR["[" *i* "]"]),
                                ")")["[" *j* "]"],
                          "]")
                  ["[" *k* "]"])

boxplot_top <- ggplot(data = boxData) +
  geom_boxplot(aes(x=ID, y=maxMaxEAR, fill = Relevant), outlier.size = 0.5) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_log10(y_label,
                labels=toxEval:::fancyNumbers,breaks=pretty_logs_new)

aop_ep <- ggplot(data = chem_sum_AOP) +
  geom_tile(aes(x=ID, y=endPoint, fill=meanEAR)) +
  theme_bw() +
  scale_x_discrete(position="top") +
  ylab("ToxCast Endpoint Name") +
  labs(fill="Mean EAR") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=7)) +
  scale_y_discrete(drop=TRUE) +
  scale_fill_gradient( guide = "legend",
                       trans = 'log',limits = c(1e-4,1),
                       low = "white", high = "steelblue",
                       breaks = c(1e-5,1e-4,1e-3,1e-2,1e-1,1),
                       labels = toxEval:::fancyNumbers2,
                       na.value = 'transparent') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        # legend.position = "none",
        plot.background = element_rect(fill = "transparent",colour = NA))


png("plots/aop_cow_leg.png", width = 1200, height = 1200, res = 142)
plot_grid(boxplot_top,  aop_ep, align = "v", nrow = 2, rel_heights = c(4/10, 6/10))
dev.off()

# How many AOPs are included, and how many are yes and maybe for relevance
endpoints <- unique(boxData$ID)
relevanceAOPs <- boxData %>% #filter(grepl("Yes",Relevant,ignore.case = TRUE)) %>%
  group_by(ID,Relevant) %>%
  summarize(medianEAR = median(maxMaxEAR)) %>%
  arrange(Relevant,desc(medianEAR)) %>%
  left_join(relevance,by=c("ID","Relevant"))

unique(yesAOPs$ID)


maybeAOPs <- boxData %>% filter(grepl("Maybe",Relevant,ignore.case = TRUE))


# site_graph <- ggplot() +
#   geom_text(data = nSites, 
#             aes(x = ID, y="# Sites", label = as.character(sitehits)), 
#             vjust = 0.5, size = 2) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank())
# 
# aop_label_graph <- ggplot() +
#   geom_text(data = nSites, 
#             aes(x = ID, y="AOP ID", label = as.character(ID)), 
#             vjust = 0.5, size = 3.5, angle = 90) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.border = element_blank())
# 
# png("plots/aop_cow.png", width = 1200, height = 1200, res = 142)
# plot_grid(site_graph, boxplot_top, aop_label_graph, aop_ep, align = "v", nrow = 4, rel_heights = c(1/20, 4/20, 1/20, 7/10))
# dev.off()


# file_name <- "landuse.csv"
# landuse <- fread(file_name)
# 
# landuse <- landuse %>%
#   mutate(site = paste0("USGS-",zeroPad(site,padTo = 8))) %>%
#   select(site, LandUse=V14)
# AOP_EAR_thres <- 0.1
# chemical_thres <- 0.01
# 1. Sum individual samples by AOP_id
# 2. Filter out EAR AOP summations by threshold EAR = 0.1
# 3. Determine percent of EAR contributed by each chemical for each AOP
# 4. Filter out chemicals with less than x% contribution (x = 1%)
# 5. List chemicals for urban and ag sites that are remaining for each AOP,
# prioritize these by % contribution to AOP EARs
# 6. Determine how many sites each of the AOPs are relevant to. Keep a list of
# sites for each AOP
# 7. At each priority site for each priority AOP, determine which chemicals are
# highest priority for each chemical from step #5


# #1:
# chemicalSummary_AOP <- chemicalSummary %>%
#   left_join(AOP, by="endPoint")
# 
# graphData <- chemicalSummary_AOP %>%
#   group_by(site, date, chnm, ID) %>%
#   summarise(sumEAR_AOP = sum(EAR)) %>%
#   filter(!is.na(ID)) 
# #2: 
# graphData <- filter(graphData, 
#                     sumEAR_AOP >= AOP_EAR_thres)
# 
# pretty_logs_new <- toxEval:::prettyLogs(graphData$sumEAR_AOP)
# y_label <- "Sum of EAR"
# 
# order_ID <- graphData %>%
#   group_by(ID) %>%
#   summarize(median = median(sumEAR_AOP)) %>%
#   arrange(median)
# 
# graphData$ID <- factor(graphData$ID, levels = order_ID$ID)
# 
# aopPlot <- ggplot(data = graphData)+
#   scale_y_log10(y_label,labels=toxEval:::fancyNumbers,breaks=pretty_logs_new) +
#   theme_bw() +
#   xlab("AOP ID") +
#   theme(plot.background = element_rect(fill = "transparent",colour = NA),
#         axis.text.y = element_text(color = "black", vjust = 0.2), 
#         axis.text.x = element_text(color = "black", vjust = 0),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0, margin = margin(-0.5,0,0,0))) +
#   geom_boxplot(aes(x=ID, y=sumEAR_AOP),lwd=0.1,outlier.size=1, fill = "steelblue")       +
#   coord_flip()
# 
# plot_info <- ggplot_build(aopPlot)
# layout_stuff <- plot_info$layout
# 
# xmin <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[1]))
# xmax <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[2]))
# ymax <- suppressWarnings(layout_stuff$panel_ranges[[1]]$y.range[2])
# 
# countNonZero <- graphData %>%
#   group_by(ID) %>%
#   summarise(nonZero = as.character(length(unique(site[sumEAR_AOP>0]))),
#             nChems = as.character(length(unique(chnm[sumEAR_AOP>0])))) %>%
#   data.frame()
# 
# aopPlot_w_labels <- aopPlot + 
#   geom_text(data=countNonZero, aes(x=ID, y=xmin,label=nonZero),size=3) +
#   geom_text(data=countNonZero, aes(x=ID, y=xmax,label=nChems),size=3) +
#   geom_text(data=data.frame(x = c(Inf,Inf), y=c(xmin,xmax), label = c("# Sites","# Chems"), stringsAsFactors = FALSE), 
#             aes(x = x,  y=y, label = label),
#             size=3)
# 
# aopPlot_w_labels
# 
# gb <- ggplot2::ggplot_build(aopPlot_w_labels)
# gt <- ggplot2::ggplot_gtable(gb)
# 
# gt$layout$clip[gt$layout$name=="panel"] <- "off"
# 
# dir.create(file.path("plots"), showWarnings = FALSE)
# png("plots/fig4_AOPid_boxplots.png", width = 800, height = 1200, res = 142)
# grid::grid.draw(gt)
# dev.off()
# 
# #3: I think this means to ignore sites?
# aop_summary <- graphData %>%
#   group_by(site, date, ID) %>%
#   summarize(total = sum(sumEAR_AOP))
# 
# fractions <- graphData %>%
#   left_join(aop_summary, by=c("site","date","ID")) %>%
#   group_by(site, date, chnm, ID) %>%
#   summarize(fraction = sumEAR_AOP/total) %>%
#   arrange(desc(ID),date) %>%
#   left_join(landuse, by="site")

###################################


# Need:
# Endpoint | AOP_vector | maxEAR | medianEAR
# AOP | endpoint_vector | maxEAR | medianEAR

# AOP_ears <- AOP %>%
#   group_by(ID) %>%
#   summarise(ep_vector = list(unique(endPoint))) 
# 
# EP_AOPs <- AOP %>%
#   group_by(endPoint) %>%
#   summarise(AOP_vector = list(unique(ID))) 
# 
# find_same_vectors <- function(test_list){
#   similar_vectors <- list()
#   for(i in 1:length(test_list)){
#     v_i <- test_list[[i]]
#     same_vectors <- sapply(test_list[-i], function(x) all(x %in% v_i))
#     if(any(same_vectors)){
#       similar_vectors[[i]] <- c(names(test_list[i]),names(same_vectors[same_vectors]))
#     }
#     
#   }
#   similar_vectors[sapply(similar_vectors, is.null)] <- NULL
#   similar_vectors <- similar_vectors[!duplicated(sapply(similar_vectors, sort))]
#   return(similar_vectors)
# }
# 
# AOPs_with_common_eps <- find_same_vectors(AOP_ears$ep_vector)
# eps_with_common_AOPs <- find_same_vectors(EP_AOPs$AOP_vector)

