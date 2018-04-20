library(toxEval)
library(dplyr)

file_name <- "OWC_data_WQ_concentrations.xlsx"
full_path <- file.path(file_name)

tox_list_concentrations <- create_toxEval(full_path)
chemicalSummary_conc <- get_chemical_summary(tox_list_concentrations)

#Trim some names:
levels(chemicalSummary_conc$Class)[levels(chemicalSummary_conc$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary_conc$Class)[levels(chemicalSummary_conc$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary_conc$Class)[levels(chemicalSummary_conc$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

graphData_conc <- graph_chem_data(chemicalSummary_conc)
graphData_conc$guide_side <- "Concentration\n[μg/L]"
