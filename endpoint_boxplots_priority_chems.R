#Plot endpoint boxplots for each of the priority chems from the manuscript
#One chemical per page of boxplots. 
#Used for analysis and formulation of the text in the discussion.

library(toxEval)
#NOTE: Add path to path_to_file!!!
#path_to_file <- 'C:/Users/srcorsi/Documents/R/win-library/3.4/toxEval/extdata/OWC_data_fromSup.xlsx' 
path_to_file <- "D:/SRCLData/Git/toxManuscript/OWC_data_fromSup.xlsx"
tox_list <- create_toxEval(path_to_file)
chem_info <- as.data.frame(tox_list$chem_info)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)


CASnums <- c("134-62-3","80-05-7","58-08-2","1912-24-9","51218-45-2","50-32-8","78-51-3","126-73-8","84852-15-3","119-61-9")
chemNames <- c(chem_info[which(chem_info$CAS %in% CASnums),2],"Fluoranthene")
chemNames[grep("pyrene",chemNames)] <- "Benzo(a)pyrene"
chemNames[grep("DEET",chemNames)] <- "DEET"
chemNames[grep("Nonylphenol",chemNames)] <- "4-Nonylphenol, branched"

filenm <- "EP_boxplots_priority_chems.pdf"
pdf(filenm, width = 8.5, height = 11)
for(i in 1:length(chemNames)){
  ######################################
  chem_name <- chemNames[i]
  ep_plot <- plot_tox_endpoints(chemical_summary, 
                                category = 'Chemical',
                                mean_logic = FALSE,
                                hit_threshold = NA,
                                title = paste('Summing EARs for',chem_name,'of a sample, taking the max of each site'),
                                filterBy = chem_name)
  print(ep_plot)
}
dev.off()
shell.exec(filenm)
