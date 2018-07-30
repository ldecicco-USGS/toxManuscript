# SI 4
library(readr)
dataDir <- "D:/LADData/RCode/toxEval_Archive/INVITRODB_V2_LEVEL5"
setwd(dataDir)
files <- list.files()

x <- read_csv(files[1], 
              col_types = list(stkc = col_double()))

x <- select(x, chnm, casn, aenm, logc_min, logc_max, modl_acc,
                   modl, actp, modl_ga, flags, hitc,gsid_rep) 

for(i in files[c(-1)]){
  subX <- read_csv(i, col_types = list(stkc = col_double())) 
  
  subFiltered <- select(subX, chnm, casn, aenm, logc_min, logc_max, modl_acc,
                        modl, actp, modl_ga, flags, hitc,gsid_rep) 
  
  x <- bind_rows(x, subFiltered)

}

# setwd("D:/LADData/RCode/toxEval_Archive/Scripts for Paper")
# 
# saveRDS(x, file="all_tox.rds")

file_name <- "D:/LADData/RCode/toxEval_Archive/Scripts for Paper/OWC_data_fromSup.xlsx"
full_path <- file.path(file_name)
tox_list <- create_toxEval(full_path)
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)
cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)
chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

total_counts <- filter(x, casn %in% tox_list$chem_info$CAS)

totals <- total_counts %>%
  filter(!is.na(modl_acc)) %>%
  group_by(casn) %>%
  summarize(Total = n(),
            min_ACC = min(modl_acc, na.rm = TRUE),
            Active = sum(hitc == 1),
            Considered = sum(hitc == 1 &
                             aenm %in% filtered_ep$endPoint ))

flag_totals <- chemicalSummary %>%
  select(CAS, endPoint) %>%
  distinct() %>%
  group_by(CAS) %>%
  summarize(Filtered = n())

totals_final <- select(tox_list$chem_info, `OWC Class` = Class, `Compound Name` = `Chemical Name`, CAS) %>%
  left_join(totals, by=c("CAS"="casn")) %>%
  left_join(flag_totals, by="CAS") %>%
  left_join(select(tox_chemicals, CAS=Substance_CASRN, mlwt = Structure_MolWt), by="CAS") %>%
  mutate(min_ACC = mlwt * (10^(min_ACC))) %>%
  arrange(`OWC Class`, desc(`Compound Name`)) 

totals_final <- totals_final[c(names(totals_final)[1:4],
                               names(totals_final)[6:8],
                               "min_ACC")]


write.csv(totals_final, file = "D:/LADData/RCode/toxEval_Archive/Scripts for Paper/tables/SI4.csv", row.names = FALSE, na = "-")
