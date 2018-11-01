################################
library(toxEval)
library(dplyr)
library(ggplot2)

file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(file_name)

tox_list <- create_toxEval(full_path)
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

sitesOrdered <- c("StLouis","Nemadji","WhiteWI","Bad","Montreal",
                  "PresqueIsle","Ontonagon","Sturgeon","Tahquamenon","Burns",
                  "IndianaHC","StJoseph","PawPaw","Kalamazoo","GrandMI",
                  "Milwaukee","Muskegon","WhiteMI","PereMarquette","Manitowoc",
                  "Manistee","Fox","Oconto","Peshtigo","Menominee",
                  "Indian","Cheboygan","Ford","Escanaba","Manistique",
                  "ThunderBay","AuSable","Rifle","Saginaw","BlackMI",
                  "Clinton","Rouge","HuronMI","Raisin","Maumee",
                  "Portage","Sandusky","HuronOH","Vermilion","BlackOH",
                  "Rocky","Cuyahoga","GrandOH","Cattaraugus","Tonawanda",
                  "Genesee","Oswego","BlackNY","Oswegatchie","Grass",
                  "Raquette","StRegis")
tox_list$chem_site$`Short Name` <- factor(tox_list$chem_site$`Short Name`,
                                          levels = sitesOrdered)
lakes_ordered <- c("Lake Superior",
                   "Lake Michigan",
                   "Lake Huron",
                   "Lake Erie",
                   "Lake Ontario")

tox_list$chem_site$site_grouping <- factor(tox_list$chem_site$site_grouping,
                                           levels=lakes_ordered)
chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)

levels(chemicalSummary$chnm)[levels(chemicalSummary$chnm) == "TDCPP"] <- "Tris(1,3-dichloro-2-propyl)phosphate"

#Trim some names:
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Antimicrobial Disinfectants"] <- "Antimicrobial"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Detergent Metabolites"] <- "Detergent"
levels(chemicalSummary$Class)[levels(chemicalSummary$Class) == "Flavors and Fragrances"] <- "Flavor/Fragrance"
