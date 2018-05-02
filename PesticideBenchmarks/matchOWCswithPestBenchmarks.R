library(dplyr)
library(toxEval)

# source("data_setup.R")
 source("data_setup_wq_benchmarks.R")

options(scipen = 10)

#Read excel file with data and benchmarks
file_name <- "OWC_data_WQ_benchmarks.xlsx"
full_path <- file.path(file_name)

tox_list <- create_toxEval(full_path)

chemicals <- as.data.frame(tox_list$chem_info)
names(chemicals)[2] <- "Chemical"


current_benchmarks <- tox_list$benchmarks
current_benchmarks <- filter(current_benchmarks,!is.na(ACC_value))
names(current_benchmarks)[4] <- "Value"
names(current_benchmarks)[2] <- "Chemical"


file_name <- "EPAAqLifeBenchmarksPest.csv"
full_path <- file.path("PesticideBenchmarks",file_name)

pest_benchmarks <- read.csv(full_path, stringsAsFactors = FALSE)

OWC_pest_benchmarks <- pest_benchmarks[which(pest_benchmarks$CASN %in% chemicals$CAS),]
OWC_pest_benchmarks_orig <- OWC_pest_benchmarks
#add column with chemical name consistent with our study
row.names(chemicals) <- chemicals$CAS
OWC_pest_benchmarks$Chemical <- chemicals[OWC_pest_benchmarks$CASN,"Chemical"]

updated_benchmarks <- tox_list$benchmarks
names(updated_benchmarks)[2] <- "Chemical"
names(updated_benchmarks)[4] <- "Value"

unique(updated_benchmarks$endPoint)
endpoints <- names(OWC_pest_benchmarks)[4:9]
groups <- c("Acute","Chronic","Acute","Chronic","Acute","Acute")


for (i in 1:length(endpoints)){
  toxTemp <- OWC_pest_benchmarks
  names(toxTemp)[3] <- "CAS"
  toxTemp$endPoint <- endpoints[i]
  toxTemp$groupCol <- groups[i]
  
  names(toxTemp)[which(names(toxTemp)==endpoints[i])] <- "Value"
  
  toxTemp <- toxTemp[,c(names(updated_benchmarks))]
  
  updated_benchmarks <- rbind(updated_benchmarks,toxTemp)
  
}

write.csv(updated_benchmarks,file="udated_benchmarks.csv",row.names = FALSE)
  
