#Explore numbers of endpoints vs active and vs min_ACC

df <- read.csv("./Tables/SI4_2_nodash.csv",stringsAsFactors = FALSE)

plot(df$Total~df$min_ACC,log="x")
plot(df$Total~df$Active,log="xy")
plot(df$Total~df$Active,log="")
plot(df$Total~df$Filtered,log="")

