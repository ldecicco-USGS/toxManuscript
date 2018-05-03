library(toxEval)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(dataRetrieval)

####################################
source(file = "data_setup.R")

AOP_EAR_thres <- 0.01
chemical_thres <- 0.01

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

file_name <- "landuse.csv"
landuse <- fread(file_name)

landuse <- landuse %>%
  mutate(site = paste0("USGS-",zeroPad(site,padTo = 8))) %>%
  select(site, LandUse=V14)

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


#1:
chemicalSummary_AOP <- chemicalSummary %>%
  left_join(AOP, by="endPoint")

graphData <- chemicalSummary_AOP %>%
  group_by(site, date, chnm, ID) %>%
  summarise(sumEAR_AOP = sum(EAR)) %>%
  filter(!is.na(ID)) 
#2: 
graphData <- filter(graphData, 
                    sumEAR_AOP >= AOP_EAR_thres)

pretty_logs_new <- toxEval:::prettyLogs(graphData$sumEAR_AOP)
y_label <- "Sum of EAR"

order_ID <- graphData %>%
  group_by(ID) %>%
  summarize(median = median(sumEAR_AOP)) %>%
  arrange(median)

graphData$ID <- factor(graphData$ID, levels = order_ID$ID)

aopPlot <- ggplot(data = graphData)+
  scale_y_log10(y_label,labels=toxEval:::fancyNumbers,breaks=pretty_logs_new) +
  theme_bw() +
  xlab("AOP ID") +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.y = element_text(color = "black", vjust = 0.2), 
        axis.text.x = element_text(color = "black", vjust = 0),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0, margin = margin(-0.5,0,0,0))) +
  geom_boxplot(aes(x=ID, y=sumEAR_AOP),lwd=0.1,outlier.size=1, fill = "steelblue")       +
  coord_flip()

plot_info <- ggplot_build(aopPlot)
layout_stuff <- plot_info$layout

xmin <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[1]))
xmax <- suppressWarnings(10^(layout_stuff$panel_ranges[[1]]$x.range[2]))
ymax <- suppressWarnings(layout_stuff$panel_ranges[[1]]$y.range[2])

countNonZero <- graphData %>%
  group_by(ID) %>%
  summarise(nonZero = as.character(length(unique(site[sumEAR_AOP>0]))),
            nChems = as.character(length(unique(chnm[sumEAR_AOP>0])))) %>%
  data.frame()

aopPlot_w_labels <- aopPlot + 
  geom_text(data=countNonZero, aes(x=ID, y=xmin,label=nonZero),size=3) +
  geom_text(data=countNonZero, aes(x=ID, y=xmax,label=nChems),size=3) +
  geom_text(data=data.frame(x = c(Inf,Inf), y=c(xmin,xmax), label = c("# Sites","# Chems"), stringsAsFactors = FALSE), 
            aes(x = x,  y=y, label = label),
            size=3)

aopPlot_w_labels

gb <- ggplot2::ggplot_build(aopPlot_w_labels)
gt <- ggplot2::ggplot_gtable(gb)

gt$layout$clip[gt$layout$name=="panel"] <- "off"

dir.create(file.path("plots"), showWarnings = FALSE)
png("plots/fig4_AOPid_boxplots.png", width = 800, height = 1200, res = 142)
grid::grid.draw(gt)
dev.off()

#3: I think this means to ignore sites?
aop_summary <- graphData %>%
  group_by(site, date, ID) %>%
  summarize(total = sum(sumEAR_AOP))

fractions <- graphData %>%
  left_join(aop_summary, by=c("site","date","ID")) %>%
  group_by(site, date, chnm, ID) %>%
  summarize(fraction = sumEAR_AOP/total) %>%
  arrange(desc(ID),date) %>%
  left_join(landuse, by="site")


