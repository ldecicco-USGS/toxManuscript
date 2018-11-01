library(toxEval)
library(dplyr)

file_name <- "OWC_data_WQ_benchmarks.xlsx"
full_path <- file.path(file_name)

tox_list_wq_benchmarks <- create_toxEval(full_path)

x <- tox_list_wq_benchmarks$benchmarks

name_fix <- tox_chemicals %>%
  select(CAS=Substance_CASRN, chnm = Substance_Name) %>%
  filter(CAS %in% unique(x$CAS))

x <- x %>%
  rename(orig_name=chnm) %>%
  left_join(name_fix, by="CAS")

x$chnm[is.na(x$chnm)] <- x$orig_name[is.na(x$chnm)]

tox_list_wq_benchmarks$benchmarks <- x

chemicalSummary_bench <- get_chemical_summary(tox_list_wq_benchmarks)
levels(chemicalSummary_bench$chnm)[levels(chemicalSummary_bench$chnm) == "TDCPP"] <- "Tris(1,3-dichloro-2-propyl)phosphate"

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

