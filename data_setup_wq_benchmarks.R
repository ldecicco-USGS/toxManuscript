library(toxEval)
library(dplyr)

file_name <- "OWC_data_WQ_benchmarks.xlsx"
full_path <- file.path(file_name)

tox_list_wq_benchmarks <- create_toxEval(full_path)

chemicalSummary_bench <- get_chemical_summary(tox_list_wq_benchmarks)

#Trim some names:
levels(chemicalSummary_bench$Class)[levels(chemicalSummary_bench$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary_bench$Class)[levels(chemicalSummary_bench$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary_bench$Class)[levels(chemicalSummary_bench$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

chemicalSummary_wqp <- filter(chemicalSummary_bench, Bio_category != "EEQ")
chemicalSummary_eeq <- filter(chemicalSummary_bench, Bio_category == "EEQ")

chemicalSummary_wqp$chnm <- factor(chemicalSummary_wqp$chnm, 
                                   levels = levels(chemicalSummary_wqp$chnm)[levels(chemicalSummary_wqp$chnm) %in% unique(as.character(chemicalSummary_wqp$chnm))])

chemicalSummary_eeq$chnm <- factor(chemicalSummary_eeq$chnm, 
                                   levels = levels(chemicalSummary_eeq$chnm)[levels(chemicalSummary_eeq$chnm) %in% unique(as.character(chemicalSummary_eeq$chnm))])

