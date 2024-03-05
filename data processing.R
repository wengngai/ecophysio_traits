library(reshape2)
library(stringr)

CNR <- read.csv("./raw_data/CN ratio.csv", header=T, stringsAsFactors = T)
Hmax <- read.csv("./raw_data/Hmax.csv", header=T, stringsAsFactors = T)
Varea <- read.csv("./raw_data/twig diameter and area.csv", header=T, stringsAsFactors = T)
vessels <- read.csv("./raw_data/vessel counts.csv", header=T, stringsAsFactors = T)
WD <- read.csv("./raw_data/wooddensity.csv", header=T, stringsAsFactors = T)
stomata <- read.csv("./raw_data/stomata.csv", header=T, stringsAsFactors = T)
lvein <- read.csv("./raw_data/leaf vein.csv", header=T, stringsAsFactors = T)
llayers <- read.csv("./raw_data/leaf layers.csv", header=T, stringsAsFactors = T)
lsoft <- read.csv("./raw_data/leaf soft traits.csv", header=T, stringsAsFactors = T)
recruitment <- read.csv("./raw_data/5cm recruitment (29 spp).csv")

### CNR, HMAX ###
# Already summarized at species mean level
summary(CNR)
summary(Hmax)

### VESSEL AREAS & DIAMETERS ###
# need to be summarized at twig then individual then species level
summary(Varea <- na.omit(Varea[1:6]))

## Area
Varea_twig <- aggregate(x = Varea["Area"], by = list(Species=Varea$Species, Individual=Varea$Individual, Twig=Varea$Twig), FUN = sum)
TS_twig <- aggregate(x = vessels$work.area, by = list(Species=vessels$Species, Individual=vessels$Individual, Twig=vessels$Twig), FUN = mean)
Varea_twig$TS_area <- TS_twig$x[match(Varea_twig$Twig, TS_twig$Twig)]
Varea_twig$VA <- log(with(Varea_twig, Area/TS_area))
Varea_indiv <- aggregate(x = Varea_twig["VA"], by = list(Species=Varea_twig$Species, Individual=Varea_twig$Individual), FUN = mean, na.rm = T)
Varea_sp <- aggregate(x = Varea_indiv["VA"], by = list(Species=Varea_indiv$Species), FUN = mean, na.rm = T)

## DH
DH <- function(x) (sum(x^4)/length(x))^(1/4)
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
vessels$T_VD <- vessels$total_vessels / vessels$work.area
vessels$VGI <- vessels$total_vessels / vessels$grouped_vessels
vessels$SV <- vessels$solitary / vessels$total_vessels
vessels$CV <- vessels$cluster / vessels$total_vessels
vessels$TV <- vessels$tangential / vessels$total_vessels
vessels$RV <- vessels$radial / vessels$total_vessels
vessels$OV <- vessels$occluded / vessels$total_vessels

vessels_indiv <- aggregate(x = vessels[12:18], by = list(Species=vessels$Species, Individual=vessels$Individual), FUN = mean, na.rm=T)
vessels_sp <- aggregate(x = vessels_indiv[3:9], by = list(Species=vessels_indiv$Species), FUN = mean, na.rm=T)

### WD ###
# need to be summarized at species level
summary(WD <- droplevels(na.omit(WD)))
WD_sp <- aggregate(x = WD[,"Density"], by = list(Species=WD$Species), FUN = mean, na.rm=T)
names(WD_sp)[2] <- "WD"

### LEAF SOFT TRAITS ###
summary(lsoft)
lsoft$species <- toupper(lsoft$species)
lsoft$SLA[which(lsoft$SLA == "#DIV/0!" | lsoft$SLA == 0)] <- NA
lsoft$ldmc[which(lsoft$ldmc == "#DIV/0!" | lsoft$ldmc == 0)] <- NA
lsoft$thickness[which(lsoft$thickness == "#DIV/0!" | lsoft$thickness == 0)] <- NA
lsoft[c("SLA", "ldmc", "thickness")] <- apply(lsoft[c("SLA", "ldmc", "thickness")], 2, as.numeric)

# need to be summarized by twig then individual then species level
lsoft_indiv <- aggregate(x = lsoft[c("SLA", "ldmc", "thickness")], by = list(Species = lsoft$species, Individual = lsoft$indiv),
          FUN = mean, na.rm = T)
lsoft_sp <- aggregate(x = lsoft_indiv[c("SLA", "ldmc", "thickness")], by = list(Species = lsoft_indiv$Species),
                         FUN = mean, na.rm = T)
names(lsoft_sp) <- c("Species", "SLA", "LDMC", "Th")

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
lvein$L_VD <- lvein$TL/lvein$Area

# need to be summarized by twig then individual then species level
lvein$Twig <- substr(lvein$Leaf, 0, 8)
lvein$Individual <- substr(lvein$Leaf, 0, 5)
lvein$Species <- substr(lvein$Leaf, 0, 3)

lvein_twig <- aggregate(x = lvein$L_VD, by = list(Species=lvein$Species, Individual=lvein$Individual, Twig=lvein$Twig), FUN = mean, na.rm=T)
lvein_indiv <- aggregate(x = lvein_twig$x, by = list(Species=lvein_twig$Species, Individual=lvein_twig$Individual), FUN = mean, na.rm=T)
lvein_sp <- aggregate(x = lvein_indiv$x, by = list(Species=lvein_indiv$Species), FUN = mean, na.rm=T)
names(lvein_sp)[2] <- "L_VD"

### LEAF LAYERS ###
summary(llayers)
llayers <- droplevels(llayers[-which(llayers$Data.collected.by==""),])
llayers$LE2.2[which(llayers$LE2.2 == "#VALUE!")] <- NA
llayers$LE2.2 <- as.numeric(as.character(llayers$LE2.2))
llayers$UE2.2[which(llayers$UE2.2 == "#VALUE!")] <- NA
llayers$UE2.2 <- as.numeric(as.character(llayers$UE2.2))

# need to be summarized by leaf then twig then individual then species level
# take out just the data columns first
llayers.dat <- llayers[,2:(length(llayers)-1)]

# by leaf
by.leaf <- as.factor(substr(names(llayers.dat), 0, 3))
llayers_leaf <- matrix(ncol=length(levels(by.leaf)), nrow=nrow(llayers.dat))
for(i in 1:nrow(llayers.dat)){
    llayers_leaf[i,] <- tapply(as.numeric(as.vector(llayers.dat[i,])), by.leaf, mean, na.rm=T)
}
colnames(llayers_leaf) <- levels(by.leaf)
llayers_leaf <- data.frame(llayers_leaf)

leaf.no <- rep(1:3, 6)
Th_123 <- cbind(
    apply(llayers_leaf[which(leaf.no == 1)], 1, sum),
    apply(llayers_leaf[which(leaf.no == 2)], 1, sum),
    apply(llayers_leaf[which(leaf.no == 3)], 1, sum)
)
Th_all <- cbind(Th_123, Th_123, Th_123, Th_123, Th_123, Th_123)
dim(Th_all); dim(llayers_leaf)
llayers_leaf <- llayers_leaf / Th_all
# outliers (data entry errors)
llayers_leaf[which(llayers_leaf[,1] > 0.1), 1] <- NA
llayers_leaf[which(llayers_leaf[,6] > 0.2), 6] <- NA
llayers_leaf[which(llayers_leaf[,14] > 0.2), 14] <- NA

# by twig
llayers_leaf$twig <- substr(llayers$Leaf, 0, 8)
llayers_twig <- melt(llayers_leaf, id.vars = "twig", measure.vars = 2:(length(llayers_leaf)-1),
     value.name = "thickness")
unique(llayers_twig$variable)
llayers_twig$variable <- as.factor(substr(as.character(llayers_twig$variable), 0, 2))
llayers_twig$Individual <- substr(as.character(llayers_twig$twig), 0, 5)
llayers_twig$Species <- substr(as.character(llayers_twig$twig), 0, 3)

# by individual then species level
llayers_indiv <- dcast(formula = Species + Individual ~ variable, 
                       value.var = "thickness", fun.aggregate = mean, na.rm = T, data = llayers_twig)
llayers_sp <- aggregate(llayers_indiv[3:8], by = list(Species = llayers_indiv$Species), FUN = mean, na.rm = T)
# check proportions seem corre:
apply((llayers_sp)[2:7], 2, mean)
# add "Th_" in front of variable names for naming
names(llayers_sp)[2:7] <- paste0("Th_", names(llayers_sp)[2:7])


### DEMOGRAPHIC RATE PARAMETERS ###
AG_parms <- read.csv("./raw_data/adult growth mod parms Sep21.csv", row.names=1)
S_parms <- read.csv("./raw_data/survival model parms Sep21.csv", row.names=1)
# rename r1 and p1 to just r and p
names(S_parms) <- c("K", "r", "p")

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
SSI <- data.frame(
    SSI = SSI$ssi.ba[match(rownames(AG_parms), rownames(SSI))],
    Species = AG_parms$sp
)

### JOINING DATA ###
traits_sp <- merge(lvein_sp, stomata_sp, by="Species")
traits_sp <- data.frame(traits_sp,
                   vessels_sp[match(traits_sp$Species, vessels_sp$Species), 2:8],
                   WD = WD_sp[match(traits_sp$Species, WD_sp$Species), "WD"],
                   lsoft_sp[match(traits_sp$Species, lsoft_sp$Species), 2:4],
                   DH = DH_sp[match(traits_sp$Species, DH_sp$Species), "DH"],
                   VA = Varea_sp[match(traits_sp$Species, Varea_sp$Species), "VA"],
                   llayers_sp[match(traits_sp$Species, llayers_sp$Species), 2:7],
                   CNR = CNR[match(traits_sp$Species, CNR$Sp), "CN.Ratio.Fresh"],
                   Hmax = Hmax[match(traits_sp$Species, Hmax$Sp), "Value"],
                   AG_parms[match(traits_sp$Species, AG_parms$sp), 1:3],
                   S_parms[match(traits_sp$Species, S_parms$sp), 1:3],
                   rec = recruitment$recruitment[match(traits_sp$Species, abbrev(recruitment$species))],
                   SSI = SSI$SSI[match(traits_sp$Species, SSI$Species)]
                   )

# define function first (at the end of script)
pairs.cor(traits_sp[2:11])
pairs.cor(traits_sp[12:24])
pairs.cor(traits_sp[25:32])

# create summary table
summary(traits_sp)
traits_sp$SD <- traits_sp$SD*1000000 # um-2 --> mm-2
traits_sp$T_VD <- traits_sp$T_VD*1000000 # um-2 --> mm-2
ranges <- data.frame(
    trait = names(traits_sp[2:24]),
    range = paste0(signif(apply(traits_sp[2:24],2,min),3), "-", signif(apply(traits_sp[2:24],2,max),3))
)
#write.csv(ranges, "./raw_data/trait ranges.csv")

### Transforming variables for PCA ###

logit <- function(p){
    if(sum(p==0)>0){
        q <- p+0.001
        return(log( q / (1 - q) ))
    } else { 
        return(log( p / (1 - p) ))
    }
}
logtrans <- function(x){
    if(sum(x==0)>0){
        offset <- min(x[which(x>0)]) / 2
        return( log(x + offset) )
    } else { 
        return(log(x))
    }
}

traits_sp[c("SV", "CV", "TV", "RV", "OV")] <-
    apply(traits_sp[c("SV", "CV", "TV", "RV", "OV")], 2, logit)
traits_sp[c("Th_LC", "Th_PM", "Th_SM", "Th_UC", "Th_UE")] <-
    apply(traits_sp[c("Th_LC", "Th_PM", "Th_SM", "Th_UC", "Th_UE")], 2, logit)
traits_sp$K <- logit(traits_sp$K)
traits_sp$p1 <- logtrans(traits_sp$p1)
traits_sp$b <- log(traits_sp$b)
traits_sp$rec <- sqrt(traits_sp$rec)
summary(traits_sp)
#write.csv(traits_sp, "combined traits_sp level_Oct23.csv")


##################################
# INVESTIGATING TRAIT PLASTICITY #
##################################

indiv_traits <- read.csv("./raw_data/pilot data for intraspecific analysis.csv")
summary(indiv_traits)
indiv_traits$indiv <- toupper(indiv_traits$indiv)
indiv_traits <- aggregate(x = indiv_traits[,7:10], by = list(Individual = indiv_traits$indiv, Hydrology = indiv_traits$hydrology), FUN = mean, na.rm = T)


# these are the data with individual level data
head(Varea_indiv)
head(llayers_indiv)
head(vessels_indiv)
head(lvein_indiv)
head(WD)
WD$individual <- sapply(str_split(WD$Twig, "-"), "[[", 1)
WD_indiv <- aggregate(x = WD["Density"], by = list(Species=WD$Species, Individual=WD$individual), FUN = mean)
head(DH_indiv)
head(stomata_indiv)

indiv_traits$VA <- Varea_indiv$VA[match(indiv_traits$Individual, Varea_indiv$Individual)]
indiv_traits$Th_SM <- llayers_indiv$SM[match(indiv_traits$Individual, llayers_indiv$Individual)]
indiv_traits$VGI <- vessels_indiv$VGI[match(indiv_traits$Individual, vessels_indiv$Individual)]
indiv_traits$T_VD <- vessels_indiv$T_VD[match(indiv_traits$Individual, vessels_indiv$Individual)]
indiv_traits$L_VD <- lvein_indiv$x[match(indiv_traits$Individual, lvein_indiv$Individual)]
indiv_traits$WD <- WD_indiv$Density[match(indiv_traits$Individual, WD_indiv$Individual)]
indiv_traits$DH <- DH_indiv$mean.diameter[match(indiv_traits$Individual, DH_indiv$Individual)]
indiv_traits$GCL <- stomata_indiv$GCL[match(indiv_traits$Individual, stomata_indiv$Individual)]
indiv_traits$SD <- stomata_indiv$SD[match(indiv_traits$Individual, stomata_indiv$Individual)]
indiv_traits$Species <- str_sub(indiv_traits$Individual, 1, 3)

library(nlme)
summary(indiv_traits)
pairs.cor(indiv_traits[,3:15])

indiv_traits$sla[which(indiv_traits$sla==0)] <- NA
indiv_traits$dbh[which(indiv_traits$dbh==0)] <- NA

hist(indiv_traits$thickness <- log(indiv_traits$thickness))
hist(indiv_traits$sla <- log(indiv_traits$sla))
hist(indiv_traits$dbh <- log(indiv_traits$dbh))
hist(indiv_traits$L_VD <- log(indiv_traits$L_VD))
hist(indiv_traits$SD <- log(indiv_traits$SD))
hist(indiv_traits$VGI <- log(indiv_traits$VGI))
hist(indiv_traits$DH <- log(indiv_traits$DH))

indiv_traits[,3:14] <- apply(indiv_traits[,3:14], 2, scale)
indiv_traits$Species <- substr(indiv_traits$Individual, 0, 3)
modlist <- list()

modlist[[1]] <- lme(sla ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[2]] <- lme(ldmc ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[3]] <- lme(thickness ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[4]] <- lme(Th_SM ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[5]] <- lme(L_VD ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[6]] <- lme(SD ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[7]] <- lme(GCL ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[8]] <- lme(WD ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[9]] <- lme(T_VD ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[10]] <- lme(VGI ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)
modlist[[11]] <- lme(DH ~ Hydrology, random = ~1|Species, data = indiv_traits, na.action = na.omit)

plot(modlist[[1]])
summary(modlist[[1]])$tTable

hydro.coef <- data.frame(coef = rep(NA, length(modlist)), upp = NA, low = NA)
#dbh.coef <- data.frame(coef = rep(NA, length(modlist)), upp = NA, low = NA)

for(i in 1:length(modlist)){
    hydro.coef$coef[i] <- summary(modlist[[i]])$tTable[2, 1]
    hydro.coef$upp[i] <- summary(modlist[[i]])$tTable[2, 1] + 1.96*summary(modlist[[i]])$tTable[2, 2]
    hydro.coef$low[i] <- summary(modlist[[i]])$tTable[2, 1] - 1.96*summary(modlist[[i]])$tTable[2, 2]
    #dbh.coef$coef[i] <- summary(modlist[[i]])$tTable[3, 1]
    #dbh.coef$upp[i] <- summary(modlist[[i]])$tTable[3, 1] + 1.96*summary(modlist[[i]])$tTable[3, 2]
    #dbh.coef$low[i] <- summary(modlist[[i]])$tTable[3, 1] - 1.96*summary(modlist[[i]])$tTable[3, 2]
}

hydro.cols <- ifelse(hydro.coef$upp < 0 | hydro.coef$low > 0, "black", "grey")
#dbh.cols <- ifelse(dbh.coef$upp < 0 | dbh.coef$low > 0, "black", "grey")

#jpeg("./outputs/trait plasticity analysis v2.jpg", width = 4, height = 7, units = "in", res = 300)
par(mar = c(4,5,2,2))
plot(seq(1, 11, 1) ~ hydro.coef$coef, pch = 16, xlim = c(-1,0.5), cex = 2, ylim = c(0.5, 12), col = hydro.cols,
     ylab = "", xlab = "Standardized effect size", yaxt = "n")
#points(seq(1.1, 11.1, 1) ~ dbh.coef$coef, pch = 17, cex = 2, ylim = c(0, 12), col = dbh.cols)
for(i in 1:length(modlist)) arrows(hydro.coef$upp[i], i, hydro.coef$low[i], i, length = 0, lwd = 2, col = hydro.cols[i])
abline(v = 0, lty = 2)
axis(side = 2, at = 1:11, labels = c("SLA", "LDMC", "Th", "Th_SM", "L_VD", "SD", "GCL", "WD", "T_VD", "VGI", "DH"), las = 1)
#legend('topleft', pch = c(16,17), legend = c("Waterlogging", "DBH"), bty = "n")
dev.off()
apply(!is.na(indiv_traits[,c(5,4,3,8,11,15,14,12,10,9,13)]), 2, sum)
# just state: btwn 79 to 109 individuals

library(MuMIn)

Rm <- function(x) round(r.squaredGLMM(x)[1], 3)
Rc <- function(x) round(r.squaredGLMM(x)[2], 3)
hydrop <- function(x) round(summary(x)$tTable[2,5], 3)
#dbhp <- function(x) round(summary(x)$tTable[3,5], 3)

Tab <- data.frame(
    trait = c("SLA", "LDMC", "Th", "Th_SM", "L_VD", "SD", "GCL", "WD", "T_VD", "VGI", "DHI"),
    hydro.p = unlist(lapply(modlist, hydrop)),
    #dbh.p = unlist(lapply(modlist, dbhp)),
    R2m = unlist(lapply(modlist, Rm)),
    R2c = unlist(lapply(modlist, Rc))
)
Tab$n <- apply(!is.na(indiv_traits[,c(5,4,3,8,11,15,14,12,10,9,13)]), 2, sum)    
#write.csv(Tab, "./outputs/Table S2.csv")

unique(indiv_traits$Species)
traits_sp$Species 
`%ni%` <- Negate(`%in%`)
traits_sp$Species[which(unique(traits_sp$Species) %ni% unique(indiv_traits$Species))]


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

