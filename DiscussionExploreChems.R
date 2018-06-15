library(toxEval)
library(readxl)

source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore DEET occurrences and range

#sum EARs by chemical per sample and take max sample per site
DEET_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS,chnm,site,shortName,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(CAS == "134-62-3")
  

boxplot(DEET_EARs$EARmax,log="y")
abline(h=10^-3)

#How many sites > 10^-3
sum(DEET_EARs$EARmax > 10^-3)

DEET_raw_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  filter(CAS == "134-62-3")

unique(DEET_raw_EARs$endPoint)

par(mar=c(12,6,3,1))
boxplot(DEET_raw_EARs$EAR~DEET_raw_EARs$endPoint,log="y",las=2)

# Deet concentrations:
DEET_conc <- concentrations %>%
  filter(Value > 0) %>%
  filter(CAS == "134-62-3" ) %>%
  group_by(SiteID) %>%
  summarize(ConcMax = max(Value))
         
boxplot(DEET_conc$ConcMax,log="y")  
range(DEET_conc$ConcMax)
summary(DEET_conc$ConcMax)

sum(DEET_conc$ConcMax>0.06)



group_by(CAS,SiteID,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) 
