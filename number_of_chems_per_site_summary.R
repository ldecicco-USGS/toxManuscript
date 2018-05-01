


EARgt0.001 <- filter(graphData,nChem > 0)

table(EARgt0.001$site_grouping)

GT15chems <- filter(graphData,nChem>=15)
