source("data_setup.R")

library(toxEval)

endPointInfo <- endPointInfo
endpointList <- c("CLD_CYP1A1_6hr","CLD_CYP1A2_6hr")
endpointListInfo <- filter(endPointInfo,assay_component_name %in% endpointList)
write.csv(endpointListInfo, file="./endPointInfo/CLD_CYP1A1_6hr.csv")
