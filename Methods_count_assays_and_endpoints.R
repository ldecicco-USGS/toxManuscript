library(toxEval)
library(readxl)

####################################
source(file = "data_setup.R")

AllEndpoints <- unique(chemicalSummary$endPoint) #255 endpoints


totalAssays <- unique(endPointInfo$assay_name)


assay_battery <- unique(substr(chemicalSummary$endPoint,1,3)) #10 assay sources

#Number of chemicals present, num in ToxCast, num with hitcalls
file_name <- "OWC_data_fromSup.xlsx"
full_path <- file.path(file_name)

df <- read_xlsx(full_path,sheet = "Data")
CASnums <- unique(df$CAS)

CASToxCast <- unique(ACC$casn)

sum(CASnums %in% CASToxCast) #54 chems in ToxCast

numChems <- length(unique(chemicalSummary$chnm)) #48 with hitcalls

endpointsPerChem <- group_by(chemicalSummary,chnm) %>%
  summarize(numEndpoints = n_distinct(endPoint))
range(endpointsPerChem$numEndpoints)



#Verifying that there are no assay sources that have gain and loss of the same assay (there are not--verified)
#Checkin on gain and loss endpoints
endPointInfo <- as.data.frame(endPointInfo)
endPointInfo2 <- endPointInfo[!(endPointInfo$assay_source_name == 
                                  "ATG" & endPointInfo$signal_direction == "loss"), ]
endPointInfo2 <- endPointInfo2[!(endPointInfo2$assay_source_name == 
                                   "NVS" & endPointInfo2$signal_direction == "gain"), ]

endPointInfo_filtered <- endPointInfo2[endPointInfo2$assay_component_endpoint_name %in% AllEndpoints,]
atgInfo <- filter(endPointInfo_filtered,assay_source_name == "ATG")

atgInfo$signal_direction

aprInfo <- filter(endPointInfo_filtered,assay_source_name == "APR")
aprEndpoints <- aprInfo$assay_component_endpoint_name
ACCApr <- ACClong[ACClong$endPoint %in% aprEndpoints,]
ACCendpointNameNo_UpDown <- sub("_up","",ACCApr$endPoint)
ACCApr$ACCendpointNameNo_UpDown <- sub("_dn","",ACCendpointNameNo_UpDown)

test <- group_by(ACCApr,chnm,ACCendpointNameNo_UpDown) %>%
  summarize(num = n_distinct(endPoint))


for (i in 1:length(assay_battery)){
  
  assayInfo <- filter(endPointInfo_filtered,assay_source_name == assay_battery[i])
  Endpoints <- assayInfo$assay_component_endpoint_name
  ACCAssay <- ACClong[ACClong$endPoint %in% Endpoints,]
  ACCendpointNameNo_UpDown <- sub("_up","",ACCAssay$endPoint)
  ACCAssay$ACCendpointNameNo_UpDown <- sub("_dn","",ACCendpointNameNo_UpDown)
  
  assayResult <- group_by(ACCAssay,chnm,ACCendpointNameNo_UpDown) %>%
    summarize(num = n_distinct(endPoint))
  
  if(i == 1) AllAssaysResult <- assayResult
  else AllAssaysResult <- rbind(assayResult,AllAssaysResult)
}

AllAssaysResult$num
