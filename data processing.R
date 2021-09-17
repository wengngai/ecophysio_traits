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
vessels$VGI <- vessels$total_vessels / vessels$grouped_vessels
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
# remove data entry error
lvein <- lvein[-which(lvein$TL > 200000),]
lvein$VD_leaf <- lvein$TL/lvein$Area

# need to be summarized by twig then individual then species level
lvein$Twig <- substr(lvein$Leaf, 0, 8)
lvein$Individual <- substr(lvein$Leaf, 0, 5)
lvein$Species <- substr(lvein$Leaf, 0, 3)

lvein_twig <- aggregate(x = lvein$VD_leaf, by = list(Species=lvein$Species, Individual=lvein$Individual, Twig=lvein$Twig), FUN = mean, na.rm=T)
lvein_indiv <- aggregate(x = lvein_twig$x, by = list(Species=lvein_twig$Species, Individual=lvein_twig$Individual), FUN = mean, na.rm=T)
lvein_sp <- aggregate(x = lvein_indiv$x, by = list(Species=lvein_indiv$Species), FUN = mean, na.rm=T)
names(lvein_sp)[2] <- "VD_leaf"





### DEMOGRAPHIC RATE PARAMETERS ###
AG_parms <- read.csv("./raw_data/adult growth mod parms Sep21.csv", row.names=1)
S_parms <- read.csv("./raw_data/surv mod parms Sep21.csv", row.names=1)

abbrev <- function(x){
    y <- gsub("\\.", " ", x)
    toupper(paste0(
        substr(lapply(strsplit(y, " "), "[[", 1), 0, 1),
        substr(sapply(strsplit(y, " "), "[[", 2), 0, 2)
    ))}
AG_parms$sp <- abbrev(rownames(AG_parms))
S_parms$sp <- abbrev(rownames(S_parms))

zeide_w_transform <- function(a, b, c, dbh) (dbh ^ b) * exp(a - dbh*c)
newlogdbh <- seq(1,50,len=200)
plot(rep(0.5, 200) ~ newlogdbh, ylim=c(0,1), xlab="DBH (log-transformed)", ylab="AGR", type="n")
library(viridis)
for(i in 1:nrow(AG_parms)){
    lines(zeide_w_transform(a = AG_parms[i,"a"], b = AG_parms[i, "b"], c = AG_parms[i, "c"], dbh = newlogdbh) ~
              newlogdbh, col = viridis(nrow(AG_parms))[i])
}

needham_w_transform <- function(K, r, p, dbh) K / (1 + exp(-r * (dbh - p) ))
newdbh <- seq(0.001, 1, len=200)
plot(rep(0.5, 200) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival probability", type="n")
for(i in 1:nrow(S_parms)){
    lines(needham_w_transform(K = S_parms[i,"K"], r = S_parms[i, "r1"], p = S_parms[i, "p1"], dbh = newdbh) ~
              newdbh, col = viridis(nrow(S_parms))[i])
}

### SSI ###
SSI <- read.csv("./raw_data/SSI Jan21.csv", row.names = 1)
SSI$sp <- abbrev(rownames(SSI))


### JOINING DATA ###
traits_sp <- merge(lvein_sp, stomata_sp, by="Species")
traits_sp <- data.frame(traits_sp,
                   vessels_sp[match(traits_sp$Species, vessels_sp$Species), 2:8],
                   WD = WD_sp[match(traits_sp$Species, WD_sp$Species), "WD"],
                   DH = DH_sp[match(traits_sp$Species, DH_sp$Species), "DH"],
                   vessel_area = Varea_sp[match(traits_sp$Species, Varea_sp$Species), "Area"],
                   CN_ratio = CNR[match(traits_sp$Species, CNR$Sp), "CN.Ratio.Fresh"],
                   Hmax = Hmax[match(traits_sp$Species, Hmax$Sp), "Value"],
                   AG_parms[match(traits_sp$Species, AG_parms$sp), 1:3],
                   S_parms[match(traits_sp$Species, S_parms$sp), 1:3],
                   SSI = SSI$ssi.ba[match(traits_sp$Species,SSI$sp)]
                   )

# define function first (at the end of script)
pairs.cor(traits_sp[2:length(traits_sp)])

### Transforming variables for PCA ###

logit <- function(p) log( p / (1 - p) )
logtrans <- function(x){
    if(sum(x==0)>0){
        offset <- min(x[which(x>0)]) / 2
        return( log(x + offset) )
    } else { 
        return(log(x))
    }
}
traits_sp$K <- logit(traits_sp$K)
traits_sp[c("prop_solitary", "prop_cluster", "prop_tangential", "prop_radial", "prop_occluded")] <-
    apply(traits_sp[c("prop_solitary", "prop_cluster", "prop_tangential", "prop_radial", "prop_occluded")], 2, logtrans)
traits_sp$p1 <- logtrans(traits_sp$p1)

pairs.cor(traits_sp[2:length(traits_sp)])

#write.csv(traits_sp, "combined traits_sp level_Sep21.csv")



## Defining pairs.cor() function
pairs.cor <- function (x,y,smooth=TRUE, digits=2,  ...)
{
    panel.cor <- function(x, y, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r.obj = cor.test(x, y,use="pairwise",...)
        r = as.numeric(r.obj$estimate)
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(txt)
        cex <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex=cex*abs(r))
    }
    panel.hist <- function(x)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col="grey")
    }
    pairs(x,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.smooth, ...)
}