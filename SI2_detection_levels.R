library(toxEval)
library(dplyr)
library(tidyr)
library(stringi)

source(file = "data_setup.R")

# Substitute max LDL or MDL for actual values:

tox_list$chem_data <- tox_list$chem_data %>%
  left_join(select(tox_list$chem_info,
                   CAS,
                   MDL = `Maximum method detection level`,
                   LDL = `Maximum laboratory reporting level`),
            by="CAS") %>%
  rowwise() %>%
  mutate(Value = max(MDL, LDL, na.rm = TRUE)) %>%
  select(SiteID, `Sample Date`, CAS, Value) %>%
  distinct()

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"

# Basically...need to swap endpoint and site to take advantage of plot_tox_boxplots code

chemicalSummary$site <- chemicalSummary$endPoint

plot_DL <- plot_tox_boxplots(chemicalSummary, "Chemical",title = "EAR per endPoint (# on left is number of affected endpoints) ")

dir.create(file.path("plots"), showWarnings = FALSE)
ggsave(plot_DL, filename = "plots/SI2_detection_levels.png", width = 5, height = 5)

