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
source(file = "MakeTitles.R")

ear_thresh <- 0.001
siteThres <- 10
# ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE)
relevance$Relevant <- MakeTitles(relevance$Relevant)

relevance <- relevance %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.) 

endpoints_sites_hits <- filter(chemicalSummary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= ear_thresh) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThres)

priority_endpoints <- endpoints_sites_hits$endPoint

boxData_max <- chemicalSummary %>%
  left_join(AOP, by="endPoint") %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)])

boxData <- boxData_max %>%
  group_by(ID,site,date) %>%
  summarize(total = sum(maxEAR))  %>%
  group_by(ID, site) %>%
  summarize(maxMaxEAR = max(total, na.rm = TRUE),
            date_picked = date[which.max(total)]) %>%
  ungroup() %>%
  filter(!is.na(ID)) 

priority_AOPs <- boxData %>%
  filter(maxMaxEAR > ear_thresh) %>%
  group_by(ID) %>%
  summarise(siteDet = n_distinct(site)) %>%
  filter(siteDet >= siteThres)

boxData_max <- filter(boxData_max, ID %in% priority_AOPs$ID)
boxData <- filter(boxData, ID %in% priority_AOPs$ID)

relevance <- relevance %>%
  filter(ID %in% priority_AOPs$ID) 

boxData <- boxData %>%
  left_join(select(relevance, ID, Relevant), by="ID") %>%
  mutate(ID = factor(ID)) 

boxData$Relevant <- factor(boxData$Relevant, levels = c("Yes","No","Maybe"))

# fractions <- boxData_tots %>%
#   left_join(boxData_max, by=c("ID","site","date")) %>%
#   right_join(boxData, by=c("ID","site","date"="date_picked")) %>%
#   group_by(ID, site, date, chnm, endPoint_used) %>%
#   summarize(fraction = maxEAR/total) %>%
#   ungroup()
# 
# endpoints_that_contribute <- fractions %>%
#   filter(fraction > ep_percent_thres) %>%
#   select(endPoint_used) %>%
#   distinct() %>%
#   pull(endPoint_used)

chem_sum_AOP <- boxData_max %>%
  filter(endPoint_used %in% priority_endpoints) %>%
  group_by(ID, endPoint_used) %>%
  summarize(maxEAR = max(maxEAR, na.rm = TRUE),
            meanEAR = mean(maxEAR, na.rm = TRUE),
            medianEAR = median(maxEAR, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(ID)) %>%
  mutate(ID = as.factor(ID),
         endPoint = as.factor(endPoint_used)) 

nSites <- boxData %>%
  select(ID, maxMaxEAR,site) %>%
  distinct() %>%
  group_by(ID) %>%
  summarize(sitehits = sum(maxMaxEAR > ear_thresh)) %>%
  ungroup() %>%
  filter(!is.na(ID),
         ID %in% priority_AOPs$ID) %>%
  mutate(ID = as.factor(ID))

chem_sum_AOP <- filter(chem_sum_AOP, ID %in% priority_AOPs$ID)
chem_sum_AOP$endPoint <- droplevels(chem_sum_AOP$endPoint )

pretty_logs_new <- toxEval:::prettyLogs(boxData$maxMaxEAR)

y_label <- bquote("max" ~ 
                    group("[", 
                          group("(",
                                sum(" max("  ~ EAR["[" *i* "]"] ~ ")"),
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
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=c("#E69F00", "#009E73", "#F0E442","#999999")) +
  scale_y_log10(y_label,
                labels=toxEval:::fancyNumbers,
                breaks=pretty_logs_new)

aop_ep <- ggplot(data = chem_sum_AOP) +
  geom_tile(aes(x=ID, y=endPoint, fill=meanEAR)) +
  theme_bw() +
  scale_x_discrete(position="top") +
  ylab("ToxCast Endpoint Name") +
  labs(fill="Mean EAR") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),#element_text(size=7),
        legend.position = "none") +
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
        plot.background = element_rect(fill = "transparent",colour = NA))


# How many AOPs are included, and how many are yes and maybe for relevance
# endpoints <- unique(boxData$ID)
# 
# relevanceAOPs <- boxData %>% #filter(grepl("Yes",Relevant,ignore.case = TRUE)) %>%
#   group_by(ID,Relevant) %>%
#   summarize(medianEAR = median(maxMaxEAR)) %>%
#   left_join(relevance,by=c("ID","Relevant")) %>%
#   group_by(ID,Relevant,Rationale) %>%
#   summarize(medianEAR = median(medianEAR))%>%
#   arrange(Relevant,desc(medianEAR)) %>%
#   filter(grepl(c("Yes|Maybe"),Relevant,ignore.case = TRUE))
# 
# write.csv(relevanceAOPs,file="relevanceAOPs.csv",row.names=FALSE)


site_graph <- ggplot() +
  geom_text(data = nSites,
            aes(x = ID, y="# Sites", label = as.character(sitehits)),
            vjust = 0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

aop_label_graph <- ggplot() +
  geom_text(data = nSites,
            aes(x = ID, y="AOP ID", label = as.character(ID)),
            vjust = 0.5, size = 3, angle = 0) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())

legend_box <- get_legend(boxplot_top + 
                         theme(legend.position = "bottom") )
legend_aop <- get_legend(aop_ep + theme(legend.position="bottom"))

png("plots/Fig4_aop_cow.png", width = 1800, height = 1200, res = 142)
plot_grid(site_graph, boxplot_top, 
          aop_label_graph, aop_ep, 
          plot_grid(legend_box, legend_aop, ncol = 2),
          align = "v", nrow = 5, 
          rel_heights = c(1/20, 8/20, 1/20, 8/20,1/10),
          labels = c("A","","","B",""))
dev.off()


######################################################################################
#Code for exploring data to be included in manuscript text
relevantEARs <- filter(boxData, grepl("yes|maybe",Relevant,ignore.case = TRUE)) 
range(relevantEARs$maxMaxEAR) # Get max EAR

# Determine which chemicals for each AOP, range of EARs, and how many sites

  #Start with chems identified from Fig1: present at 10 or more sites with EARmaxChem > 10-3 at 5 or more sites
priority_chems <- read.csv("priority_chems.csv",stringsAsFactors = TRUE)

AOPs_by_priority_chem <- filter(boxData_max,CAS %in% priority_chems$CAS) %>%
  left_join(relevance,by=("ID")) %>%
  filter(maxEAR > 0) %>%
  filter(grepl("yes|maybe",Relevant,ignore.case = TRUE))%>%
  group_by(ID,chnm,CAS) %>%
  summarize(nSites = n_distinct(site),
            maxEAR = max(maxEAR),
            maxEndpoint = endPoint[which.max(maxEAR)]) %>%
  filter(maxEAR > 0.001) %>%
  arrange(ID,nSites)

unique(AOPs_by_priority_chem$chnm)
priority_chems[which(!priority_chems$CAS %in% AOPs_by_priority_chem$CAS),]


# Are AOPs 26, 52, and 53 significantly diff than other relevant AOPs?
relevantData <- boxData %>% 
#  filter(maxMaxEAR > 0) %>%
  filter(!Relevant == "No")
  
large <- pull(relevantData[relevantData$ID %in% c(29,52,53),],maxMaxEAR)
small <- pull(relevantData[!relevantData$ID %in% c(29,52,53),],maxMaxEAR)
plot(x=percent_rank(large), large,pch=20,col="blue",log="y")
points(x=percent_rank(small),small,col="orange")

boxplot(large, log="y")
wilcox.test(large,small)
t.test(log(large+10e-9),log(small+10^-9))
t.test(large,small)

large <- pull(relevantData[relevantData$ID %in% c(29),],maxMaxEAR)
small <- pull(relevantData[relevantData$ID %in% c(63),],maxMaxEAR)
plot(x=percent_rank(large), large,pch=20,col="blue",log="y")
points(x=percent_rank(small),small,col="orange")

wilcox.test(large,small)
t.test(large,small)

AOPsInData <- as.character(unique(relevantData$ID))
relevant_AOP_in_data_info <-relevance %>%
  filter(ID %in% AOPsInData) %>%
  group_by(ID,Relevant,Rationale) %>%
  summarize(endpoints = paste(endPoint,collapse="; ")) %>%
  inner_join(AOP_crosswalk,by=c("ID"="AOP..")) %>%
  group_by(ID,Relevant,Rationale,endpoints) #%>%
    
    AOP_crosswalk$Component.Endpoint.Name <- as.factor(AOP_crosswalk$Component.Endpoint.Name )
  relevant_AOP_in_data_info <-relevance %>%
    filter(ID %in% AOPsInData) %>%
    group_by(ID,Relevant,Rationale) %>%
    summarize(endpoints = paste(endPoint,collapse="; ")) %>%
    inner_join(AOP_crosswalk,by=c("ID"="AOP..")) %>%
    mutate(Component.Endpoint.Name =as.factor(Component.Endpoint.Name)) %>%
    group_by(ID,Relevant,Rationale,endpoints) %>%
    summarize(endpointID_crosswalk = paste(Assay.Endpoint.ID,collapse = ";"),
 #             endpoints_cwalk = paste(unique(Component.Endpoint.Name),collpse = "; "))
            AOP.Title = paste(unique(AOP.Title),collapse = "; "),
            KE. = paste(unique(KE.),collapse = "; "),
            Key.Event.Name = paste(unique(Key.Event.Name),collapse = "; "),
            KeyEvent.Type = paste(unique(KeyEvent.Type),collapse = "; "))
            
  write.csv(relevant_AOP_in_data_info,file="relevant_AOPs_in_study.csv",row.names = FALSE)

length(unique(relevant_AOP_in_data_info$ID))

relevantAOPInfo <- relevance[relevance$ID %in% as.character(unique(relevantData$ID)),] %>%
  inner_join(AOP_crosswalk,by=c("ID"="AOP.."))



# AOPs_by_chnm_all <- left_join(boxData_max,relevance,by=("ID")) %>%
#   filter(maxEAR > 0) %>%
#   filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
#   group_by(ID,chnm) %>%
#   summarize(nSites = n_distinct(site),
#             maxEAR = max(maxEAR)) %>%
#   arrange(ID,nSites)
#   
# 
# AOPs_by_chnm_thresh <- left_join(boxData_max,relevance,by=("ID")) %>%
#   filter(maxEAR > ear_thresh) %>%
#   filter(grepl("yes|maybe",Relevant,ignore.case = TRUE)) %>%
#   group_by(ID,chnm) %>%
#   summarize(nSites = n_distinct(site),
#             maxEAR = max(maxEAR)) %>%
#   filter(nSites > siteThres) %>%
#   arrange(ID,nSites)


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

