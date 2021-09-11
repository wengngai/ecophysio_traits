library(reshape2)

CNR <- read.csv("./raw_data/CN ratio.csv", header=T, stringsAsFactors = T)
Hmax <- read.csv("./raw_data/Hmax.csv", header=T, stringsAsFactors = T)
Varea <- read.csv("./raw_data/twig diameter and area.csv", header=T, stringsAsFactors = T)
vessels <- read.csv("./raw_data/vessel counts.csv", header=T, stringsAsFactors = T)
WD <- read.csv("./raw_data/wooddensity.csv", header=T, stringsAsFactors = T)

# CNR and Hmax are already summarized at species mean level
summary(CNR)
summary(Hmax)

# Vessel areas need to be summarized at twig then individual then species level
summary(Varea <- Varea[1:6])

# Vessel counts need to be summarized at individual then species level
summary(vessels <- na.omit(vessels[1:11]))
vessels_twig <- aggregate(total_vessels + grouped_vessels + solitary + cluster + tangential + radial + occluded ~ 
              Species + Twig, data = vessels, FUN = mean)
vessels_twig <- aggregate(x = vessels[5:11], by = list(vessels$Species, vessels$Individual, vessels$Twig), FUN = mean)
vessels_indiv <- aggregate(x = vessels_twig[4:10], by = list(vessels_twig$Group.1, vessels_twig$Group.2), FUN = mean)
vessels_sp <- aggregate(x = vessels_indiv[3:9], by = list(vessels_indiv$Group.1), FUN = mean)
vessels_sp

# WDs need to be summarized at species level
summary(WD <- droplevels(na.omit(WD)))
with(WD, tapply(Density, Species, mean))


