# Need to run fig 2 analysis before running the following script

#Determine numbers for writing info into text

###
#How many rivers per lake have EAR > XX


graphData %>%
  filter(nChem>0) %>%
  group_by(site_grouping) %>%
  summarize(numSites = n_distinct(`Short Name`))

graphData %>%
  group_by(site_grouping) %>%
  summarize(numSites = n_distinct(`Short Name`))

#Num sites with more than 5 chems
countThresh <- 10
graphData %>%
  filter(nChem>=countThresh) %>%
  group_by(site_grouping) %>%
  summarize(numSites = n_distinct(`Short Name`)) %>%
  summarize(numSites = sum(numSites))



# 3 with 15 or more chems
# 11 with 10 or more chems
# 21 with 5 or more chems