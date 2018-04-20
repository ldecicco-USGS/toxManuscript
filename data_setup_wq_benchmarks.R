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

graphData_wq <- graph_chem_data(chemicalSummary_wqp)
graphData_wq$guide_side <- "WQ Guidelines\nMaximum Toxicity Quotient per Site"
graphData_eeq <- graph_chem_data(chemicalSummary_eeq)
graphData_eeq$guide_side <- "Estradiol Equivalents"
