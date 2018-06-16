library(toxEval)
library(readxl)

directory <- "./ECOTOX"


## CAFFEINE ## --------------------------------------------------------------
caffeine <- read_xlsx("./ECOTOX/ECOTOXCaffeine.xlsx")
names(caffeine) <-gsub(" ","",names(caffeine))
names(caffeine) <-gsub(")","",names(caffeine))
names(caffeine) <-gsub("\\(","",names(caffeine))
names(caffeine) <-gsub("%","",names(caffeine))

caffeine$Conc1MeanStandardized <- as.numeric(caffeine$Conc1MeanStandardized)

caffeineugL <- caffeine[caffeine$Conc1UnitsStandardized == "ug/L",]
par(mar=c(16,6,3,1))

boxplot(caffeineugL$Conc1MeanStandardized~caffeineugL$SpeciesGroup,log="y",las=2,cex.axis=0.6)


# Now look at study results

source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore caffeine occurrences and range

# Caffeine concentrations:
Caffeine_conc <- concentrations %>%
  filter(Value > 0) %>%
  filter(CAS == "58-08-2" ) %>%
  group_by(SiteID) %>%
  summarize(ConcMax = max(Value))

sum(Caffeine_conc$ConcMax > 0.1)

boxplot(Caffeine_conc$ConcMax,log="y")
abline(h=10^-3)


# Caffeine EARs
source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore Caffeine occurrences and range

#sum EARs by chemical per sample and take max sample per site
Caffeine_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS,chnm,site,shortName,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(CAS == "58-08-2")


boxplot(Caffeine_EARs$EARmax,log="y")
abline(h=10^-3)

#How many sites > 10^-3
sum(Caffeine_EARs$EARmax > 10^-3)



## TBP ## --------------------------------------------------------------######################################
tbp <- read_xlsx("./ECOTOX/ECOTOXTBP.xlsx")
names(tbp) <-gsub(" ","",names(tbp))
names(tbp) <-gsub(")","",names(tbp))
names(tbp) <-gsub("\\(","",names(tbp))
names(tbp) <-gsub("%","",names(tbp))

tbp$Conc1MeanStandardized <- as.numeric(tbp$Conc1MeanStandardized)

tbpugL <- tbp[tbp$Conc1UnitsStandardized == "ug/L",]
par(mar=c(16,6,3,1))

boxplot(tbpugL$Conc1MeanStandardized~tbpugL$SpeciesGroup,log="y",las=2,cex.axis=0.6)


# Now look at study results

source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore tbp occurrences and range


# tbp concentrations:
tbp_conc <- concentrations %>%
  filter(Value > 0) %>%
  filter(CAS == "126-73-8" ) %>%
  group_by(SiteID) %>%
  summarize(ConcMax = max(Value))

sum(tbp_conc$ConcMax > 01)

boxplot(tbp_conc$ConcMax,log="y")
abline(h=10^-3)


# tbp EARs

#Explore tbp occurrences and range

#sum EARs by chemical per sample and take max sample per site
tbp_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS,chnm,site,shortName,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(CAS == "126-73-8")


boxplot(tbp_EARs$EARmax,log="y")
abline(h=10^-3)

#How many sites > 10^-3
sum(tbp_EARs$EARmax > 10^-3)

ACCTBP <- get_ACC("126-73-8")


## TBOEP ## --------------------------------------------------------------######################################
tboep <- read_xlsx("./ECOTOX/ECOTOXtboep.xlsx")
names(tboep) <-gsub(" ","",names(tboep))
names(tboep) <-gsub(")","",names(tboep))
names(tboep) <-gsub("\\(","",names(tboep))
names(tboep) <-gsub("%","",names(tboep))

tboep$Conc1MeanStandardized <- as.numeric(tboep$Conc1MeanStandardized)

tboepugL <- tboep[tboep$Conc1UnitsStandardized == "ug/L",]
par(mar=c(16,6,3,1))

boxplot(tboepugL$Conc1MeanStandardized~tboepugL$SpeciesGroup,log="y",las=2,cex.axis=0.6)


# Now look at study results

source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore tboep occurrences and range


# tboep concentrations:
tboep_conc <- concentrations %>%
  filter(Value > 0) %>%
  filter(CAS == "78-51-3" ) %>%
  group_by(SiteID) %>%
  summarize(ConcMax = max(Value))

sum(tboep_conc$ConcMax > 01)

boxplot(tboep_conc$ConcMax,log="y")
abline(h=10^-3)


# tboep EARs

#Explore tboep occurrences and range

#sum EARs by chemical per sample and take max sample per site
tboep_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS,chnm,site,shortName,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(CAS == "78-51-3")


boxplot(tboep_EARs$EARmax,log="y")
abline(h=10^-3)

#How many sites > 10^-3
sum(tboep_EARs$EARmax > 10^-3)


## IndChem ## --------------------------------------------------------------

chem <- "4-NP"
CASnum <- "84852-15-3"
IndChem <- read.delim("./ECOTOX/ECOTOX4NP.rep",sep = "|")
#IndChem <- read_xlsx("./ECOTOX/ECOTOXCaffeine.xlsx")
names(IndChem) <-gsub(" ","",names(IndChem))
names(IndChem) <-gsub(")","",names(IndChem))
names(IndChem) <-gsub("\\(","",names(IndChem))
names(IndChem) <-gsub("%","",names(IndChem))
names(IndChem) <-gsub("\\.","",names(IndChem))
names(IndChem) <-gsub("X","",names(IndChem))

IndChem$Conc1MeanStandardized <- as.numeric(IndChem$Conc1MeanStandardized)

IndChemugL <- IndChem[IndChem$Conc1UnitsStandardized == "ug/L",]
par(mar=c(16,6,3,1))

unique(IndChem$Conc1UnitsStandardized)
table(IndChem$Conc1UnitsStandardized)

boxplot(IndChemugL$Conc1MeanStandardized~IndChemugL$SpeciesGroup,log="y",las=2,cex.axis=0.6)


# Now look at study results

source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore IndChem occurrences and range

# IndChem concentrations:
IndChem_conc <- concentrations %>%
  filter(Value > 0) %>%
  filter(CAS == CASnum ) %>%
  group_by(SiteID) %>%
  summarize(ConcMax = max(Value))

sum(IndChem_conc$ConcMax > 0.1)

boxplot(IndChem_conc$ConcMax,log="y")
abline(h=10^-3)


# IndChem EARs
source("data_setup.R")
concentrations <- read_xlsx("OWC_data_fromSup.xlsx",sheet=1)
names(concentrations)[2] <- "date"

#Explore IndChem occurrences and range

#sum EARs by chemical per sample and take max sample per site
IndChem_EARs <- chemicalSummary %>%
  filter(EAR > 0) %>%
  group_by(CAS,chnm,site,shortName,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(CAS,chnm,site,shortName) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(CAS == CASnum)


boxplot(IndChem_EARs$EARmax,log="y")
abline(h=10^-3)

#How many sites > 10^-3
sum(IndChem_EARs$EARmax > 10^-3)


