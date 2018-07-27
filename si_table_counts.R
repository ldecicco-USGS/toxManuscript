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

setwd("D:/LADData/RCode/toxEval_Archive/Scripts for Paper")

saveRDS(x, file="all_tox.rds")

source(file = "data_setup.R")


total_counts <- filter(x, casn %in% tox_list$chem_info$CAS)

totals <- total_counts %>%
  group_by(casn) %>%
  summarize(starting_count = n(),
            total_with_hitc = sum(hitc == 1),
            total_after_base_filter = sum(hitc == 1 &
                                            aenm %in% filtered_ep$endPoint ))

flag_totals <- chemicalSummary %>%
  select(CAS, endPoint) %>%
  distinct() %>%
  group_by(CAS) %>%
  summarize(total_after_flag_removed = n())

totals <- select(tox_list$chem_info, `OWC Class` = Class, `Compound Name` = `Chemical Name`, CAS) %>%
  left_join(totals, by=c("CAS"="casn")) %>%
  left_join(flag_totals, by="CAS") 

totals <- totals %>%
  arrange(`OWC Class`, desc(`Compound Name`))

totals$starting_count[is.na(totals$starting_count)] <- 0

write.csv(totals, file = "tables/SI4.csv", row.names = FALSE, na = "-")
