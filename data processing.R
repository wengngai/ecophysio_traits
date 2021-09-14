library(reshape2)

CNR <- read.csv("./raw_data/CN ratio.csv", header=T, stringsAsFactors = T)
Hmax <- read.csv("./raw_data/Hmax.csv", header=T, stringsAsFactors = T)
Varea <- read.csv("./raw_data/twig diameter and area.csv", header=T, stringsAsFactors = T)
vessels <- read.csv("./raw_data/vessel counts.csv", header=T, stringsAsFactors = T)
WD <- read.csv("./raw_data/wooddensity.csv", header=T, stringsAsFactors = T)
stomata <- read.csv("./raw_data/stomata.csv", header=T, stringsAsFactors = T)
lvein <- read.csv("./raw_data/leaf vein.csv", header=T, stringsAsFactors = T)

### CNR, HMAX ###
# Already summarized at species mean level
summary(CNR)
summary(Hmax)

### VESSEL AREAS & DIAMETERS ###
# need to be summarized at twig then individual then species level
summary(Varea <- na.omit(Varea[1:6]))

## Area
Varea_twig <- aggregate(x = Varea["Area"], by = list(Species=Varea$Species, Individual=Varea$Individual, Twig=Varea$Twig), FUN = mean)
Varea_indiv <- aggregate(x = Varea_twig["Area"], by = list(Species=Varea_twig$Species, Individual=Varea_twig$Individual), FUN = mean)
Varea_sp <- aggregate(x = Varea_indiv["Area"], by = list(Species=Varea_indiv$Species), FUN = mean)

## DH
DH <- function(x) sum((x^4)/length(x))^(1/4)
Varea$mean.diameter <- apply(Varea[c("Major.diameter","Minor.diameter")], 1, mean)
DH_twig <- aggregate(x = Varea["mean.diameter"], 
                     by = list(Species=Varea$Species, Individual=Varea$Individual, Twig=Varea$Twig),
                     FUN = DH)
DH_indiv <- aggregate(x = DH_twig["mean.diameter"], by = list(Species=DH_twig$Species, Individual=DH_twig$Individual), FUN = mean)
DH_sp <- aggregate(x = DH_indiv["mean.diameter"], by = list(Species=DH_indiv$Species), FUN = mean)
names(DH_sp)[2] <- "DH" 

### VESSEL COUNTS ###
# need to be summarized at individual then species level
summary(vessels <- na.omit(vessels[1:11]))
vessels$VD_twig <- vessels$total_vessels / vessels$work.area
vessels$prop_grouped <- vessels$grouped_vessels / vessels$total_vessels
vessels$prop_solitary <- vessels$solitary / vessels$total_vessels
vessels$prop_cluster <- vessels$cluster / vessels$total_vessels
vessels$prop_tangential <- vessels$tangential / vessels$total_vessels
vessels$prop_radial <- vessels$radial / vessels$total_vessels
vessels$prop_occluded <- vessels$occluded / vessels$total_vessels

vessels_indiv <- aggregate(x = vessels[12:18], by = list(Species=vessels$Species, Individual=vessels$Individual), FUN = mean, na.rm=T)
vessels_sp <- aggregate(x = vessels_indiv[3:9], by = list(Species=vessels_indiv$Species), FUN = mean, na.rm=T)

### WD ###
# need to be summarized at species level
summary(WD <- droplevels(na.omit(WD)))
WD_sp <- aggregate(x = WD[,"Density"], by = list(Species=WD$Species), FUN = mean, na.rm=T)
names(WD_sp)[2] <- "WD"

### STOMATA ###
summary(stomata)
stomata <- droplevels(stomata[-which(stomata$Data.collected.by==""),])
stomata$SD1 <- as.numeric(as.character(stomata$SD1))
stomata$SD <- apply(stomata[,c("SD1", "SD2", "SD3", "SD4", "SD5")], 1, mean, na.rm=T)
stomata$GCL <- apply(stomata[,c("GCL1", "GCL2", "GCL3", "GCL4", "GCL5", "GCL6", "GCL7", "GCL8", "GCL9", "GCL10")], 1, mean, na.rm=T)
# remove data entry error
stomata <- stomata[-which(stomata$GCL > 100),]
# outliers which may be potentially problematic (Gironniera) (don't remove for now)
stomata <- stomata[-which(stomata$GCL > 30),]
plot(SD ~ GCL, data=stomata)

# need to be summarized by twig then individual then species level
stomata$Twig <- substr(stomata$Leaf, 0, 8)
stomata$Individual <- substr(stomata$Leaf, 0, 5)
stomata$Species <- substr(stomata$Leaf, 0, 3)

stomata_twig <- aggregate(x = stomata[c("GCL", "SD")], by = list(Species=stomata$Species, Individual=stomata$Individual, Twig=stomata$Twig), FUN = mean, na.rm=T)
stomata_indiv <- aggregate(x = stomata_twig[c("GCL", "SD")], by = list(Species=stomata_twig$Species, Individual=stomata_twig$Individual), FUN = mean, na.rm=T)
stomata_sp <- aggregate(x = stomata_indiv[c("GCL", "SD")], by = list(Species=stomata_indiv$Species), FUN = mean, na.rm=T)

### LEAF VEIN ###
summary(lvein)
lvein <- droplevels(lvein[-which(lvein$Data.collected.by==""),])
lvein <- lvein[-which(lvein$Area=="#VALUE!"),]
lvein$Area <- as.numeric(as.character(lvein$Area))
lvein$TL1 <- as.numeric(as.character(lvein$TL1))
lvein$TL <- apply(lvein[c("TL1", "TL2", "TL3")], 1, mean)
lvein$VD_leaf <- lvein$TL/lvein$Area

# need to be summarized by twig then individual then species level
lvein$Twig <- substr(lvein$Leaf, 0, 8)
lvein$Individual <- substr(lvein$Leaf, 0, 5)
lvein$Species <- substr(lvein$Leaf, 0, 3)

lvein_twig <- aggregate(x = lvein$VD_leaf, by = list(Species=lvein$Species, Individual=lvein$Individual, Twig=lvein$Twig), FUN = mean, na.rm=T)
lvein_indiv <- aggregate(x = lvein_twig$x, by = list(Species=lvein_twig$Species, Individual=lvein_twig$Individual), FUN = mean, na.rm=T)
lvein_sp <- aggregate(x = lvein_indiv$x, by = list(Species=lvein_indiv$Species), FUN = mean, na.rm=T)
names(lvein_sp)[2] <- "VD_leaf"


### JOINING DATA ###
traits_sp <- merge(lvein_sp, stomata_sp, by="Species")
traits_sp <- data.frame(traits_sp,
                   vessels_sp[match(traits_sp$Species, vessels_sp$Species), 2:8],
                   WD = WD_sp[match(traits_sp$Species, WD_sp$Species), "WD"],
                   DH = DH_sp[match(traits_sp$Species, DH_sp$Species), "DH"],
                   vessel_area = Varea_sp[match(traits_sp$Species, Varea_sp$Species), "Area"],
                   CN_ratio = CNR[match(traits_sp$Species, CNR$Sp), "CN.Ratio.Fresh"],
                   Hmax = Hmax[match(traits_sp$Species, Hmax$Sp), "Value"]
                   )
traits_sp












