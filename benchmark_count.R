library(dplyr)
bm <- read.csv("updated_benchmarks.csv",stringsAsFactors = FALSE)
#data.frame(bm$Value,as.numeric(bm$Value))
bm2 <- bm %>%
  filter(!is.na(Value)) %>%
  filter(endPoint != "EEQ")

unique(bm2$CAS)
