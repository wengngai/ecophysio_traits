library(vegan)
traits_sp <- read.csv("combined traits_sp level_Oct21.csv", row.names = 1, header = T)
summary(traits_sp)

#####################
# PHYLOGENETIC TREE #
#####################

library(brranching)
#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(phytools)

spp_list <- read.csv("./raw_data/species list.csv")
V.spp_list <- data.frame(
  species = gsub(" ", "_", spp_list$Species),
  genus = sapply(strsplit(spp_list$Species, " "), "[[", 1),
  family = spp_list$Family
)
tree <- phylo.maker(V.spp_list, scenarios = c("S1", "S2", "S3"))
plotTree(tree$scenario.3)

library(stringr)
abbrev.tip <- function(x){
    names <- str_split(x, "_")
    abbrev.phylo <- toupper(paste0(substr(sapply(names, "[[", 1), 1, 1), substr(sapply(names, "[[", 2), 1, 2)))
    abbrev.phylo[which(abbrev.phylo=="PCL")] <- "ACL"
    return(abbrev.phylo)
}

# sort traits_sp into order of phylo tree first
traits_sp <- traits_sp[match(abbrev.tip(tree$scenario.3$tip.label), traits_sp$Species),]

# EITHER OR:
tree$tip.label <- abbrev.tip(tree$tip.label)
#tree$tip.label <- spp_list$Species[match(abbrev.tip(tree$tip.label), spp_list$Sp)]
#tree$tip.label <- spp_list$long_name[match(abbrev.tip(tree$tip.label), spp_list$Sp)]

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Phylogenetic tree.pdf", width = 8, height = 8)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Phylogenetic tree.jpg", width = 8, height = 8, units = "in", res = 300)
plot(tree$scenario.3, no.margin=TRUE,
     tip.color = ifelse(traits_sp$SSI > 0.67, "#538A95", ifelse(traits_sp$SSI < 0.33, "#FC7753", "#B1A792")))
dev.off()

#######################
# DIMENSION REDUCTION #
#        (PCAs)       #
#######################

### Demographic params ###
names(traits_sp)
PCA.demog <- rda(traits_sp[25:32], scale=T)
summary(PCA.demog)$cont
rownames(PCA.demog$CA$u) <- traits_sp$Species

#pdf("./outputs/demog PCA1-2.pdf", width=5, height=5)
biplot(PCA.demog, choices=c(1,2), cex=3, cex.lab=1.5) 
dev.off()
# PC1 = fast-slow growth continuum (lower PC1 = faster growth rates) 
# PC2 = small-large stature continuum (higher PC2 = larger stature)
# invert so that higher PC1 = faster growth
PCA.demog$CA$u[,1] <- -1 * PCA.demog$CA$u[,1]
PCA.demog$CA$v[,1] <- -1 * PCA.demog$CA$v[,1]

#pdf("./outputs/demog PCA1-3.pdf", width=5, height=5)
biplot(PCA.demog, choices=c(4,3))
dev.off()
# PC3 = recruitment-growth/survival tradeoff (higher PC3 = more recruitment, higher survival)
# PC4 = ??

demog.PC1 <- PCA.demog$CA$u[,1]
demog.PC2 <- PCA.demog$CA$u[,2]
demog.PC3 <- PCA.demog$CA$u[,3]
demog.PC4 <- PCA.demog$CA$u[,4]


### Traits ###
PCA.traits <- rda(traits_sp[2:24], scale = T)
summary(PCA.traits)$cont
rownames(PCA.traits$CA$u) <- traits_sp$Species

#pdf("./outputs/traits PCA1-2.pdf", width=5, height=5)
biplot(PCA.traits, choices = c(1,2), display = "species")
dev.off()
# PC1: water acquisitive(+)-conservative(-) spectrum
# PC2: hydraulic safety versus hydraulic efficiency?
biplot(PCA.traits, choices = c(3,4), display = "species")
# PC3: sclerophylly?
# PC4: resource acquisitive(+)-conservative spectrum
biplot(PCA.traits, choices = c(4,6), display = "species")

# invert so that higher tPC4 = more acquisitive
PCA.traits$CA$u[,4] <- -1 * PCA.traits$CA$u[,4]
PCA.traits$CA$v[,4] <- -1 * PCA.traits$CA$v[,4]

# RDA of traits against SSI and demographic params
rda.ssi <- rda(traits_sp[2:24] ~ traits_sp$SSI + demog.PC1 + demog.PC2 + demog.PC3 + demog.PC4, scale = T)
anova(rda.ssi, by="margin")

### Investigating the effect of traits (PCs) on demographic and environmental niches

# extract trait PCs
tPC1 <- PCA.traits$CA$u[,1]
tPC2 <- PCA.traits$CA$u[,2]
tPC3 <- PCA.traits$CA$u[,3]
tPC4 <- PCA.traits$CA$u[,4]
tPC5 <- PCA.traits$CA$u[,5]
tPC6 <- PCA.traits$CA$u[,6]
tPC7 <- PCA.traits$CA$u[,7]
tPC8 <- PCA.traits$CA$u[,8]

# PCA of everything (demog params and traits)
PCA.all <- rda(traits_sp[2:33], scale=T)
biplot(PCA.all, choices = c(1,2), display = "species")
summary(PCA.all)$cont

######################
# PGLS: SSI ~ traits #
######################

library(MuMIn)
library(nlme)
library(ape)

logit <- function(x){
  x[which(x>0.99)] <- 0.99
  x[which(x<0.01)] <- 0.01
  return(log(x/(1-x)))
}
hist(logit(traits_sp$SSI),breaks=10)
traits_sp$logitSSI <- logit(traits_sp$SSI)

summary(max.ssi <- gls(model = logitSSI ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
                       data = traits_sp, correlation=corBrownian(phy = tree$scenario.3), method="ML"))
head(dr.ssi <- dredge(max.ssi, extra = "R^2"))
dr.ssi[which(dr.ssi$df==2),]
summary(best.dssi <- get.models(dr.ssi, subset=1)[[1]])

biplot(PCA.traits, choices=c(6,3), display = "species") 

# Calculate R2
library(rr2)
# First create null models
m0.ssi <- gls(logitSSI ~ 1, correlation=corBrownian(phy = tree$scenario.3), data = traits_sp, method="ML")
R2.lik(get.models(dr.ssi, subset=2)[[1]], m0.ssi)
R2.lik(get.models(dr.ssi, subset=3)[[1]], m0.ssi)
AICc(m0.ssi)

# AIC table
(AICtable <- data.frame(dr.ssi[dr.ssi$delta <= 2,]))
AICtable$R.2[1] <- R2.lik(best.dssi, m0.ssi) #R2 = 60.6%
AICtable$R.2[2] <- R2.lik(get.models(dr.ssi, subset=2)[[1]], m0.ssi)
AICtable$R.2[3] <- R2.lik(get.models(dr.ssi, subset=3)[[1]], m0.ssi)

#write.csv(AICtable, "D:\\Dropbox\\Functional Traits Project\\Figures\\AIC table SSI only.csv")
  
############################
# TRAIT-TRAIT CORRELATIONS #
############################

traits_datonly <- traits_sp[,2:24]
traits_datonly_scaled <- apply(traits_datonly, 2, scale)
pairwise <- data.frame(t(combn(ncol(traits_datonly_scaled),2)))
names(pairwise) <- c("y", "x")

for(i in 1:nrow(pairwise)){
  y <- traits_datonly_scaled[, pairwise[i,"y"] ]
  x <- traits_datonly_scaled[, pairwise[i,"x"] ]
  mod <- gls(y ~ x, correlation=corBrownian(phy = tree$scenario.3), method="ML")
  pairwise$y.name[i] <- names(traits_datonly)[ pairwise[i,"y"] ]
  pairwise$x.name[i] <- names(traits_datonly)[ pairwise[i,"x"] ]
  pairwise$coef[i] <- coef(mod)["x"]
  pairwise$p.value[i] <- summary(mod)$tTable["x", "p-value"]
  pairwise$AIC[i] <- summary(mod)$AIC
}
pairwise[which(pairwise$p.value < 0.05),]

library(spaa)
# have to convert to pairwise to a dist object
sig.pairs <- pairwise[which(pairwise$p.value < 0.05),]
pairwise.dist <- list2dist(data.frame(x = c(sig.pairs$x.name, sig.pairs$y.name),
                                      y = c(sig.pairs$y.name, sig.pairs$x.name),
                                      coef = c(sig.pairs$coef, sig.pairs$coef)
))
pairwise.dist[is.na(pairwise.dist)] <- 0

# qgraph plots
#library(qgraph)
#attr(pairwise.dist, "Labels")
#traitgroups <- list(
#  leaf=c(1,3,5,11:19),
#  vessels=c(2,6:10,20:23),
#  wood=c(24,4)
#)
#
#Q <- qgraph(pairwise.dist, groups = traitgroups, 
#       minimum = 0, vsize = 5, 
#       legend = T, borders = F, details = T)

library(circlize)

ordered.traitnames <- c(
  "SLA", "LDMC", "Th", "CNR",
  "Th_LC", "Th_LE", "Th_SM", "Th_PM", "Th_UE", "Th_UC",
  "GCL", "SD", "L_VD",
  "VA", "T_VD", "RV", "CV", "SV", "OV", "TV", "VGI", "DH",
  "WD"
)
groupnames <- c(rep("leaf", 13), rep("wood", 10))
grid.col <- ifelse(groupnames=="leaf", "#DBD56E", "#403D58")
names(groupnames) <- ordered.traitnames
names(grid.col) <- ordered.traitnames

sig.pairs <- pairwise[c("x.name","y.name","coef")]
sig.pairs$coef[which(pairwise$p.value > 0.05)] <- 0
links.col <- ifelse(sig.pairs$coef>0, "#66D7D1", ifelse(sig.pairs$coef==0, "#F2EFEA", "#FC7753"))
alpha.sig <- (0.05 - pairwise$p.value)/0.1
alpha.sig[which(alpha.sig < 0)] <- 0
for(i in 1:length(links.col)) links.col[i] <- adjustcolor(links.col[i], alpha.f = alpha.sig[i])

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait correlations chord.pdf", height=6, width=10)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait correlations chord.jpg", height=6, width=10, units="in", res=300)
chordDiagram(sig.pairs, symmetric = T, 
             order = ordered.traitnames, group = groupnames, annotationTrack = c("name", "grid"), 
             annotationTrackHeight = mm_h(c(3, 3)),
             grid.col = grid.col, col = links.col, transparency = 0.75
)
legend(x= 1.1, y = -0.5, title = "Correlation", bty = "n",
       legend=c("Positive", "Negative"), fill = adjustcolor(c("#66D7D1", "#FC7753"), alpha.f = 0.25))
legend(x= 1.1, y = -0.2, title = "Organ", bty = "n",
       legend=c("Leaf", "Wood"), fill = c("#DBD56E", "#403D58"))
dev.off()
circos.clear()

#################
# VISUALIZATION #
#################

### SSI correlations with trait PCs 1, 3, 6
summary(best.dssi)
inv.logit <- function(y) exp(y) / (1 + exp(y))

# predictions
newtPC1 <- seq(-0.35, 0.4, 0.01)
dpred1 <- predict(best.dssi, newdata = data.frame(tPC1 = newtPC1, tPC3 = 0, tPC6 = 0), se.fit = T)
newtPC3 <- seq(-0.5, 0.5, 0.01)
dpred3 <- predict(best.dssi, newdata = data.frame(tPC3 = newtPC3, tPC1 = 0, tPC6 = 0), se.fit = T)
newtPC6 <- seq(-0.3, 0.6, 0.01)
dpred6 <- predict(best.dssi, newdata = data.frame(tPC6 = newtPC6, tPC1 = 0, tPC3 = 0), se.fit = T)

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\SSI correlations.pdf", height=10, width=4)
par(mfrow = c(3, 1), mar = c(6,7,1,2), mgp = c(4,1,0))
# a) PC1
plot(traits_sp$SSI ~ tPC1, type="n", xlim = c(-0.35, 0.4),
     ylab = "", xlab = "Trait PC1\n(Leaf sclerophylly dimension)",
     cex.lab = 1.4)
polygon(with(dpred1, inv.logit(c(fit + se.fit, rev(fit - se.fit)))) ~ c(newtPC1,rev(newtPC1)), border = F, col = "#66D7D16e")
lines(inv.logit(dpred1$fit) ~ newtPC1, lwd = 4, col = "white")
text(traits_sp$SSI ~ tPC1, labels=traits_sp$Species, 
     cex = ifelse(traits_sp$Species %in% c("LMU", "TWA"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("LMU", "TWA"), "black", "grey60"))
mtext(side = 3, text = " a)", cex = 1.2, line = -1.5, adj = 0)

# b) PC3
plot(traits_sp$SSI ~ tPC3, type="n", xlim = c(-0.5, 0.5),
     ylab = "", xlab = "Trait PC3\n(Leaf fleshiness dimension)",
     cex.lab = 1.4)
polygon(with(dpred3, inv.logit(c(fit + se.fit, rev(fit - se.fit)))) ~ c(newtPC3,rev(newtPC3)), border = F, col = "#66D7D16e")
lines(inv.logit(dpred3$fit) ~ newtPC3, lwd = 4, col = "white")
text(traits_sp$SSI ~ tPC3, labels=traits_sp$Species, 
     cex = ifelse(traits_sp$Species %in% c("PAX", "ACL"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("PAX", "ACL"), "black", "grey60"))
mtext(side = 3, text = " b)", cex = 1.2, line = -1.5, adj = 0)

# c) PC6
plot(traits_sp$SSI ~ tPC6, type="n", xlim = c(-0.3, 0.6),
     ylab = "", xlab = "Trait PC6\n(Vein density dimension)",
     cex.lab = 1.4)
polygon(with(dpred6, inv.logit(c(fit + se.fit, rev(fit - se.fit)))) ~ c(newtPC6,rev(newtPC6)), border = F, col = "#66D7D16e")
lines(inv.logit(dpred6$fit) ~ newtPC6, lwd = 4, col = "white")
text(traits_sp$SSI ~ tPC6, labels=traits_sp$Species, 
     cex = ifelse(traits_sp$Species %in% c("MBE", "PEC"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("MBE", "PEC"), "black", "grey60"))
mtext(side = 3, text = " c)", cex = 1.2, line = -1.5, adj = 0)
mtext(side = 2, text = "SSI (Swamp association)", cex = 1.3, outer = T, line = -3)

dev.off()

### Trait effect on demog 
# read in again untransformed params
AG_parms <- read.csv("./raw_data/adult growth mod parms Sep21.csv", row.names=1)
S_parms <- read.csv("./raw_data/surv mod parms Sep21.csv", row.names=1)

AG_parms$sp <- abbrev(rownames(AG_parms))
S_parms$sp <- abbrev(rownames(S_parms))

zeide_w_transform <- function(a, b, c, dbh) (dbh ^ b) * exp(a - dbh*c)
needham_w_transform <- function(K, r, p, dbh) K / (1 + exp(-r * (dbh - p) ))

# assign colors
col.hi1 <- "#66D7D1"
col.lo1 <- "#B1A792"
col.hi2 <- "#538A95"
col.lo2 <- "#FC7753"
col.hi3 <- "#403D58"
col.lo3 <- "#F7B39F"

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Demog effects plots.pdf", height=8, width=8)
layout(matrix(c(1,2,1,3), ncol = 2))
par(mar = c(4.5,5.5,1,1))

## Trait PC2 effect on AG
# ranked: RCI, PEC, PCO (but take AFR instead) are lowest, BPA, HCR, MGI are highest
#tPC2[order(tPC2)]
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.7), xlab="DBH (cm)", 
     ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Rhodamnia cinerea","a"], b = AG_parms["Rhodamnia cinerea", "b"], c = AG_parms["Rhodamnia cinerea", "c"], dbh = newdbh) ~
        newdbh, col = col.lo1, lwd = 3)
lines(zeide_w_transform(a = AG_parms["Pternandra echinata","a"], b = AG_parms["Pternandra echinata", "b"], c = AG_parms["Pternandra echinata", "c"], dbh = newdbh) ~
        newdbh, col = col.lo1, lwd = 3)
lines(zeide_w_transform(a = AG_parms["Horsfieldia crassifolia","a"], b = AG_parms["Horsfieldia crassifolia", "b"], c = AG_parms["Horsfieldia crassifolia", "c"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 3)
lines(zeide_w_transform(a = AG_parms["Macaranga gigantea","a"], b = AG_parms["Macaranga gigantea", "b"], c = AG_parms["Macaranga gigantea", "c"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 3)
legend('topright', bty = "n", title = "Trait PC2",
       legend = c("High: HCR, MGI", "Low: RCI, PEC"),
       lwd = 3, col = c(col.hi1, col.lo1), cex = 1.2)
mtext(side = 3, adj = 0, line = -1.5, text = " a)", cex = 1.5)

## Trait PC4 effect on survival
# ranked: RCI, AAN, PPI are lowest, TWA, AFR (but take XFL isntead), ASY are highest
plot(traits_sp$rec^2 ~ tPC4, type = "n", ylab = expression(paste("BSR (stems ", year^-1, " ", plot^-1, " ", m^-2, ")")),
     xlab = "Trait PC4", cex.lab = 1.5, xlim = c(-0.4, 0.6))
text(traits_sp$rec^2 ~ tPC4, labels = traits_sp$Species, 
     col = ifelse(traits_sp$Species %in% c("RCI", "AAN"), col.lo2,
                  ifelse(traits_sp$Species %in% c("ASY", "XFL"), col.hi2, "grey")),
     cex = ifelse(traits_sp$Species %in% c("RCI", "AAN", "ASY", "XFL"), 1.4, 0.8),
     font = ifelse(traits_sp$Species %in% c("RCI", "AAN", "ASY", "XFL"), 2, 1)
)
mtext(side = 3, adj = 0, line = -1.5, text = " b)", cex = 1.5)

newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Rhodamnia cinerea","K"], r = S_parms["Rhodamnia cinerea", "r1"], p = S_parms["Rhodamnia cinerea", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo2, lwd = 3)
lines(needham_w_transform(K = S_parms["Alstonia angustifolia","K"], r = S_parms["Alstonia angustifolia", "r1"], p = S_parms["Alstonia angustifolia", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo2, lwd = 3)
lines(needham_w_transform(K = S_parms["Aporosa symplocoides","K"], r = S_parms["Aporosa symplocoides", "r1"], p = S_parms["Aporosa symplocoides", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi2, lwd = 3)
lines(needham_w_transform(K = S_parms["Xanthophyllum flavescens","K"], r = S_parms["Xanthophyllum flavescens", "r1"], p = S_parms["Xanthophyllum flavescens", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi2, lwd = 3)
legend('bottomright', bty = "n", title = "Trait PC4",
       legend = c("High: ASY, XFL", "Low: RCI, AAN"),
       lwd = 3, col = c(col.hi2, col.lo2), cex = 1.2)
mtext(side = 3, adj = 0, line = -1.5, text = " c)", cex = 1.5)
dev.off()

### Trait PCA biplots

summary(PCA.traits)$cont
screeplot(PCA.traits)

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait PCA.pdf", height=14, width=6)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait PCA.jpg", height=14, width=6, units="in", res=300)

# a) PC1 - PC2
par(mfrow = c(3,1), mar = c(5,5.5,2,2))
plot(x = PCA.traits$CA$v[,1]*1.25, y = PCA.traits$CA$v[,2]*1.25, type = "n",
     ylab = "Trait PC2 (19.5% variance explained)", xlab = "Trait PC1 (22.5% variance explained)",
     cex.lab = 1.5)
for(i in 1:nrow(PCA.traits$CA$v)){
  arrowcol <- ifelse(i %in% c(4:11,15:16), "#DBD56E", "#403D58")
  arrows(0, 0, PCA.traits$CA$v[i,1], PCA.traits$CA$v[i,2], lwd = 2, col = arrowcol)
  text(x = PCA.traits$CA$v[i,1]*1.2, y = PCA.traits$CA$v[i,2]*1.2, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.traits$CA$v)[i], cex = 1.2)
}
mtext(text = " a)", side = 3, adj = 0, cex = 1.5, line = -2)

# b) PC3 - PC4
plot(x = PCA.traits$CA$v[,3]*1.25, y = PCA.traits$CA$v[,4]*1.25, type = "n",
     ylab = "Trait PC4 (9.9% variance explained)", xlab = "Trait PC3 (13.3% variance explained)",
     cex.lab = 1.5)
for(i in 1:nrow(PCA.traits$CA$v)){
  arrowcol <- ifelse(i %in% c(4:11,15:16), "#DBD56E", "#403D58")
  arrows(0, 0, PCA.traits$CA$v[i,3], PCA.traits$CA$v[i,4], lwd = 2, col = arrowcol)
  text(x = PCA.traits$CA$v[i,3]*1.2, y = PCA.traits$CA$v[i,4]*1.2, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.traits$CA$v)[i], cex = 1.2)
}
mtext(text = " b)", side = 3, adj = 0, cex = 1.5, line = -2)

# c) PC5 - PC6
plot(x = PCA.traits$CA$v[,5]*1.25, y = PCA.traits$CA$v[,6]*1.25, type = "n",
     xlab = "Trait PC5 (8.0% variance explained)", ylab = "Trait PC6 (5.5% variance explained)",
     cex.lab = 1.5)
for(i in 1:nrow(PCA.traits$CA$v)){
  arrowcol <- ifelse(i %in% c(4:11,15:16), "#DBD56E", "#403D58")
  arrows(0, 0, PCA.traits$CA$v[i,5], PCA.traits$CA$v[i,6], lwd = 2, col = arrowcol)
  text(x = PCA.traits$CA$v[i,5]*1.2, y = PCA.traits$CA$v[i,6]*1.2, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.traits$CA$v)[i], cex = 1.2)
}
mtext(text = " c)", side = 3, adj = 0, cex = 1.5, line = -2)

legend('bottomright', legend = c("Leaf traits", "Wood traits"), lwd = 2, col = c("#403D58", "#DBD56E"), bty = "n")

dev.off()

### Trait-Demog correlations

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait-demog correlations.pdf", height=10, width=5)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait-demog correlations.jpg", height=10, width=5, units="in", res=900)
layout(matrix(c(1,2,4,1,3,4), ncol = 2))
par(mar = c(6,7,2,2), mgp = c(4,1,0))
# a) higher demog.PC1 (growth rate)  driven by lower tPC2 (greater hydraulic investment)
newtPC2 <- seq(min(tPC2), max(tPC2), 0.01)
apred <- predict(best.dPC1, newdata = data.frame(tPC2 = newtPC2, tPC4 = 0), se.fit = T)
plot(demog.PC1 ~ tPC2, type = "n", las = 1,
     ylab = "Demographic PC1\n(fast-slow growth dimension)",
     xlab = "Trait PC2\n(Hydraulic investment dimension)",
     cex.lab = 1.4)
polygon(with(apred, c(fit + se.fit, rev(fit - se.fit))) ~ c(newtPC2,rev(newtPC2)), border = F, col = "#F7B39F6e")
lines(apred$fit ~ newtPC2, lwd = 4, col = "white")
text(demog.PC1 ~ tPC2, labels = names(tPC2), 
     cex = ifelse(traits_sp$Species %in% c("RCI", "GNE"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("RCI", "GNE"), "black", "grey60"))
mtext(side = 3, text = " a)", cex = 1.2, line = -1.5, adj = 0)

# b) higher demog.PC2 (stature, survivability/longevity) enabled by higher investment in spongy mesophyl (tPC1)
newtPC1 <- seq(min(tPC1), max(tPC1), 0.01)
bpred <- predict(best.dPC2, newdata = data.frame(tPC1 = newtPC1, tPC4 = 0), se.fit = T)
plot(demog.PC2 ~ tPC1, type = "n", las = 1, xlim = c(-0.45,0.35),
     xlab = "",
     ylab = "Demographic PC2\n(small-large stature dimension)",
     cex.lab = 1.4)
polygon(with(bpred, c(fit + se.fit, rev(fit - se.fit))) ~ c(newtPC1,rev(newtPC1)), border = F, col = "#F7B39F6e")
lines(bpred$fit ~ newtPC1, lwd = 4, col = "white")
text(demog.PC2 ~ tPC1, labels = names(tPC1), 
     cex = ifelse(traits_sp$Species %in% c("SCE", "MGI"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("SCE", "MGI"), "black", "grey60"))
mtext(side = 3, text = " b)", cex = 1.2, line = -1.5, adj = 0)

# c) higher demog.PC3 (recruitment, survival) weakly driven by higher tPC1 (thinner cuticles/epidermises, thicker spongy mesophylls)
newtPC1 <- seq(min(tPC1), max(tPC1), 0.01)
cpred <- predict(best.dPC3, newdata = data.frame(tPC1 = newtPC1, tPC2 = 0, tPC4 = 0, tPC6 = 0), se.fit = T)
plot(demog.PC3 ~ tPC1, type = "n", xlim = c(-0.45,0.35),
     xlab = "",
     ylab = "Demographic PC3\n(recruitment/survival-growth  dimension)",
     cex.lab = 1.4)
polygon(with(cpred, c(fit + se.fit, rev(fit - se.fit))) ~ c(newtPC1,rev(newtPC1)), border = F, col = "#F7B39F6e")
lines(cpred$fit ~ newtPC1, lwd = 4, col = "white")
text(demog.PC3 ~ tPC1, labels = names(tPC1), 
     cex = ifelse(traits_sp$Species %in% c("SCE", "MGI"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("SCE", "MGI"), "black", "grey60"))
mtext(side = 3, text = " c)", cex = 1.2, line = -1.5, adj = 0)

mtext(side = 1, outer = T, line = -27, text ="Trait PC1\n(Spongy mesophyll dimension)", adj = 0.6)

# d) swamp adaptation driven by higher vein densities with less occlusion (higher PC6)
newtPC6 <- seq(min(tPC6), max(tPC6), 0.01)
dpred <- predict(best.dssi, newdata = data.frame(tPC6 = newtPC6, tPC4 = 0, tPC5 = 0), se.fit = T)
plot(traits_sp$SSI ~ tPC6, type="n",
     ylab = "SSI\n(Swamp association)",
     xlab = "Trait PC6\n(Vein density dimension)",
     cex.lab = 1.4)
polygon(with(dpred, c(fit + se.fit, rev(fit - se.fit))) ~ c(newtPC6,rev(newtPC6)), border = F, col = "#66D7D16e")
lines(dpred$fit ~ newtPC6, lwd = 4, col = "white")
text(traits_sp$SSI ~ tPC6, labels=traits_sp$Species, 
     cex = ifelse(traits_sp$Species %in% c("MBE", "PEC"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("MBE", "PEC"), "black", "grey60"))
mtext(side = 3, text = " d)", cex = 1.2, line = -1.5, adj = 0)

dev.off()


#########################
#      OBSOLETIZED:     #
# PGLS FOR DEMOG TRAITS #
#########################

library(MuMIn)
library(nlme)
library(ape)

### demog.PC1
summary(max.dPC1 <- gls(demog.PC1 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
                        correlation=corBrownian(phy = tree$scenario.3), method="ML"))
#summary(lm.dPC1 <- lm(demog.PC1 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, na.action = na.fail))
head(dr.dPC1 <- dredge(max.dPC1, extra = "R^2", m.lim = c(0,4)))
#(dr.dPC1.lm <- dredge(lm.dPC1, extra = "R^2", m.lim = c(0,4)))
summary(best.dPC1 <- get.models(dr.dPC1, subset=1)[[1]])

# higher demog.PC1 (growth rate) driven by higher tPC4 (resource acquisitive)
plot(demog.PC1 ~ tPC4, type = "n"); text(demog.PC1 ~ tPC4, labels = names(tPC4))

# higher demog.PC1 (growth rate)  driven by higher tPC2 (greater hydraulic investment)
plot(demog.PC1 ~ tPC2, type = "n"); text(demog.PC1 ~ tPC2, labels = names(tPC2))

#pdf("./outputs/traits PCA2-4.pdf", width=5, height=5)
biplot(PCA.traits, choices=c(4,2), display = "species") 
dev.off()

### demog.PC2
summary(max.dPC2 <- gls(demog.PC2 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
                        correlation=corBrownian(phy = tree$scenario.3), method="ML"))
#summary(lm.dPC2 <- gls(demog.PC2 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, na.action = na.fail))
head(dr.dPC2 <- dredge(max.dPC2, extra = "R^2", m.lim = c(0,4)))
#(dr.dPC2.lm <- dredge(lm.dPC2, extra = "R^2", m.lim = c(0,4)))
summary(best.dPC2 <- get.models(dr.dPC2, subset=1)[[1]])
# higher demog.PC2 (stature, survivability/longevity) enabled by higher investment in water transport. but correlations weak
plot(demog.PC2 ~ tPC6, type = "n"); text(demog.PC2 ~ tPC6, labels = names(tPC6))
plot(demog.PC2 ~ tPC7, type = "n"); text(demog.PC2 ~ tPC7, labels = names(tPC7))

biplot(PCA.traits, choices=c(6,7), display = "species") 

### demog.PC3
summary(max.dPC3 <- gls(demog.PC3 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
                        correlation=corBrownian(phy = tree$scenario.3), method="ML"))
head(dr.dPC3 <- dredge(max.dPC3, extra = "R^2", m.lim = c(0,4)))
summary(best.dPC3 <- get.models(dr.dPC3, subset=1)[[1]])
# higher demog.PC3 (recruitment, survival) weakly driven by lower tPC1 (thinner cuticles/epidermises, thicker spongy mesophylls)
plot(demog.PC3 ~ tPC1, type = "n"); text(demog.PC3 ~ tPC1, labels = names(tPC1))
# higher demog.PC3 (recruitment, survival) driven by lower tPC2 (lesser hydraulic investment)
plot(demog.PC3 ~ tPC2, type = "n"); text(demog.PC3 ~ tPC2, labels = names(tPC2))
# higher demog.PC3 (recruitment, survival) driven by higher tPC4 (more acquisitive leaf resource acquisition strategy)
plot(demog.PC3 ~ tPC4, type = "n"); text(demog.PC3 ~ tPC4, labels = names(tPC4))
# higher demog.PC3 (recruitment, survival) weakly driven by lower tPC3 
# (thicker, more sclerophyllous leaves with lower palisade meso and LDMC) but this trait only found in top model
plot(demog.PC3 ~ tPC3, type = "n"); text(demog.PC3 ~ tPC3, labels = names(tPC3))

biplot(PCA.traits, choices=c(1,2), display = "species")
biplot(PCA.traits, choices=c(3,4), display = "species")

biplot(PCA.demog, choices=c(3,4))

### demog.PC4
summary(max.dPC4 <- gls(demog.PC4 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
                        correlation=corBrownian(phy = tree$scenario.3), method="ML"))
(dr.dPC4 <- dredge(max.dPC4, extra = "R^2"))

# swamp adaptation weakly driven by lower PC4 (more conservative leaf resource acquisition strategy)
plot(traits_sp$SSI ~ tPC4, type="n"); text(traits_sp$SSI ~ tPC4, labels=traits_sp$Species)
# swamp adaptation enabled by leaf sclerophylly (lower PC3)
plot(traits_sp$SSI ~ tPC3, type="n"); text(traits_sp$SSI ~ tPC3, labels=traits_sp$Species)
# swamp adaptation driven by higher vein densities with less occlusion (lower PC6)
plot(traits_sp$SSI ~ tPC6, type="n"); text(traits_sp$SSI ~ tPC6, labels=traits_sp$Species)

# are demog.PCs and SSI correlated?
summary(ssi.mod <- gls(model = SSI ~ demog.PC1 + demog.PC2 + demog.PC3 + demog.PC4,
                       data = traits_sp, correlation=corBrownian(phy = tree$scenario.3), method="ML"))
(dr.ssi.demo <- dredge(ssi.mod))
summary(get.models(dr.ssi.demo, subset = 1)[[1]])
# swamp species tend to be more large statured, long-lived species (maybe just the result of plot selection)
plot(traits_sp$SSI ~ demog.PC2, type="n"); text(traits_sp$SSI ~ demog.PC2, labels=traits_sp$Species)

# AIC table
AICtable <- rbind(
  rep(NA, length(dr.dPC1)),
  data.frame(dr.dPC1[dr.dPC1$delta < 2,]),
  rep(NA, length(dr.dPC1)),
  data.frame(dr.dPC2[dr.dPC2$delta < 2,]),
  rep(NA, length(dr.dPC1)),
  data.frame(dr.dPC3[dr.dPC3$delta < 2,]),
  rep(NA, length(dr.dPC1)),
  data.frame(dr.ssi[dr.ssi$delta < 2,]))

# Calculate R2
library(rr2)
# First create null models
m0.dPC1 <- gls(demog.PC1 ~ 1, correlation=corBrownian(phy = tree$scenario.3), method="ML")
m0.dPC2 <- gls(demog.PC2 ~ 1, correlation=corBrownian(phy = tree$scenario.3), method="ML")
m0.dPC3 <- gls(demog.PC3 ~ 1, correlation=corBrownian(phy = tree$scenario.3), method="ML")
m0.ssi <- gls(SSI ~ 1, correlation=corBrownian(phy = tree$scenario.3), data = traits_sp, method="ML")

# demog PC1
n1 <- nrow(dr.dPC1[dr.dPC1$delta < 2,])
for(i in 1:n1){
  mod <- get.models(dr.dPC1, subset = i)[[1]]
  AICtable$R.2[i+1] <- R2.lik(mod, m0.dPC1)
}
# demog PC2
n2 <- nrow(dr.dPC2[dr.dPC2$delta < 2,])
for(i in 1:n2){
  mod <- get.models(dr.dPC2, subset = i)[[1]]
  AICtable$R.2[i+n1+2] <- R2.lik(mod, m0.dPC2)
}
# demog PC3
n3 <- nrow(dr.dPC3[dr.dPC3$delta < 2,])
for(i in 1:n3){
  mod <- get.models(dr.dPC3, subset = i)[[1]]
  AICtable$R.2[i+n1+n2+3] <- R2.lik(mod, m0.dPC3)
}

# SSI
AICtable$R.2[16] <- R2.lik(best.dssi, m0.ssi)

AICtable

#write.csv(AICtable, "D:\\Dropbox\\Functional Traits Project\\Figures\\AIC table.csv")

### OBSOLETIZED: 
### Demographic PCA plot ###

# assign colors
col.hi1 <- "#66D7D1"
col.lo1 <- "#B1A792"
col.hi2 <- "#538A95"
col.lo2 <- "#FC7753"
col.hi3 <- "#403D58"
col.lo3 <- "#F7B39F"

summary(PCA.demog)$cont

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Demog PCA.pdf", height=12, width=9)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Demog PCA.jpg", height=12, width=9, units="in", res=300)
layout(matrix(c(1,1,2,2,3,1,1,2,2,4,5,6,7,8,9), ncol = 3))
par(mar = c(4.5,5.5,1,1))
# a) PCA of PC2 ~ PC1
plot(demog.PC2 ~ demog.PC1, type = "n", ylim = c(-0.45,0.45), xlim = c(-0.6, 0.4),
     ylab = "Demographic PC2 (21.8% vriance explained)", xlab = "Demographic PC1 (33.5% variance explained)", cex.lab = 1.5)
for(i in 1:8){
  if(i==1|i==8) arrowcol <- "#DBD56E" else if(i<=4) arrowcol <- "#403D58" else arrowcol <- "#8E8963" 
  arrows(0, 0, PCA.demog$CA$v[i,1]/2, PCA.demog$CA$v[i,2]/2, lwd = 2, col = arrowcol)
  text(x = PCA.demog$CA$v[i,1]/1.8, y = PCA.demog$CA$v[i,2]/1.8, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.demog$CA$v)[i], cex = 1.2 )
}
text(demog.PC2 ~ demog.PC1, labels = traits_sp$Species,
     col = ifelse(traits_sp$Species %in% c("LMU","SCE"), col.hi2, 
                  ifelse(traits_sp$Species %in% c("ASY","EMA"), col.lo2,
                         ifelse(traits_sp$Species %in% c("PPO", "GNE"), col.hi1,
                                ifelse(traits_sp$Species %in% c("RCI", "AAN"), col.lo1, "grey")))),
     cex = ifelse(traits_sp$Species %in% c("LMU", "ASY", "SCE", "EMA", "PPO", "GNE", "RCI", "AAN"), 1.2, 0.7)
)
legend("bottomleft", bty = "n", lwd = 2, col=c("#403D58","#8E8963", "#DBD56E"),
       legend = c("Growth model parameters: a, b, c (panels c & h)", 
                  "Survival model parameters: K, r, p (panels d, f & i)", 
                  "Other parameters: Hmax (panel e), rec (panel g)")
)
mtext(side = 3, line = -1.3, text = " a)", adj = 0)

# b) PCA of PC3 ~ PC1
plot(demog.PC3 ~ demog.PC1, type = "n", ylim = c(-0.45,0.45), xlim = c(-0.6, 0.4),
     ylab = "Demographic PC3 (13.3% variance explained)", xlab = "Demographic PC1 (33.5% variance explained)", cex.lab = 1.5)
for(i in 1:8){
  if(i==1|i==8) arrowcol <- "#DBD56E" else if(i<=4) arrowcol <- "#403D58" else arrowcol <- "#8E8963" 
  arrows(0, 0, PCA.demog$CA$v[i,1]/2, PCA.demog$CA$v[i,3]/2, lwd = 2, col = arrowcol)
  text(x = PCA.demog$CA$v[i,1]/1.8, y = PCA.demog$CA$v[i,3]/1.8, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.demog$CA$v)[i], cex = 1.2 )
}
text(demog.PC3 ~ demog.PC1, labels = traits_sp$Species,
     col = ifelse(traits_sp$Species %in% c("AFR","GAX"), col.hi3, 
                  ifelse(traits_sp$Species %in% c("ACL","MGI"), col.lo3, "grey")),
     cex = ifelse(traits_sp$Species %in% c("AFR", "GAX", "ACL", "MGI"), 1.2, 0.7)
)
#legend("bottomleft", bty = "n", lwd = 2, col=c("#403D58","#8E8963", "#DBD56E"),
#       legend = c("Growth model parameters: a, b, c (panel g)", 
#                  "Survival model parameters: K, r, p (panel h)", 
#                  "Recruitment: rec")
#)
mtext(side = 3, line = -1.3, text = " b)", adj = 0)

# c) Growth (Demog PC1)
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.7), xlab="DBH (cm)", ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Prunus polystachya","a"], b = AG_parms["Prunus polystachya", "b"], c = AG_parms["Prunus polystachya", "c"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Gironniera nervosa","a"], b = AG_parms["Gironniera nervosa", "b"], c = AG_parms["Gironniera nervosa", "c"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Rhodamnia cinerea","a"], b = AG_parms["Rhodamnia cinerea", "b"], c = AG_parms["Rhodamnia cinerea", "c"], dbh = newdbh) ~
        newdbh, col = col.lo1, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Alstonia angustifolia","a"], b = AG_parms["Alstonia angustifolia", "b"], c = AG_parms["Alstonia angustifolia", "c"], dbh = newdbh) ~ 
        newdbh, col = col.lo1, lwd = 2)
#text(x = rep(20,4), y = c(0.65, 0.48, 0.26, 0.17), labels = c("PPO", "GNE", "AAN", "RCI"), col=c(col.hi1, col.hi1, col.lo1, col.lo1))
mtext(side = 3, line = -1.3, text = " c)", adj = 0)

# d) Survival (Demog PC1)
newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Gironniera nervosa","K"], r = S_parms["Gironniera nervosa", "r1"], p = S_parms["Gironniera nervosa", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 2)
lines(needham_w_transform(K = S_parms["Prunus polystachya","K"], r = S_parms["Prunus polystachya", "r1"], p = S_parms["Prunus polystachya", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi1, lwd = 2)
lines(needham_w_transform(K = S_parms["Alstonia angustifolia","K"], r = S_parms["Alstonia angustifolia", "r1"], p = S_parms["Alstonia angustifolia", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo1, lwd = 2)
lines(needham_w_transform(K = S_parms["Rhodamnia cinerea","K"], r = S_parms["Rhodamnia cinerea", "r1"], p = S_parms["Rhodamnia cinerea", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo1, lwd = 2)
legend('bottomright', bty = "n", legend = c("High PC1: GNE, PPO", "Low PC1: AAN, RCI"), col = c(col.hi1, col.lo1), lwd = 2)
mtext(side = 3, line = -1.3, text = " d)", adj = 0)

# e) Hmax (Demog PC2)
par(mar = c(4.5,7,1,1))
h.bar <- traits_sp$Hmax[match(c("LMU", "SCE", "EMA", "ASY"),traits_sp$Species)]
barplot(h.bar, beside = T, col = c(col.hi2, col.hi2, col.lo2, col.lo2),
        ylab = "Hmax (m)", xlab = "Species", ylim = c(0,50),
        names.arg = c("LMU", "SCE", "EMA", "ASY"), cex.lab = 1.5, cex.names = 1.2)
box()
mtext(side = 3, line = -1.3, text = " e)", adj = 0)

# f) Survival (Demog PC2)
newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Lophopetalum multinervium","K"], r = S_parms["Lophopetalum multinervium", "r1"], p = S_parms["Lophopetalum multinervium", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi2, lwd = 2)
lines(needham_w_transform(K = S_parms["Strombosia ceylanica","K"], r = S_parms["Strombosia ceylanica", "r1"], p = S_parms["Strombosia ceylanica", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi2, lwd = 2)
lines(needham_w_transform(K = S_parms["Aporosa symplocoides","K"], r = S_parms["Aporosa symplocoides", "r1"], p = S_parms["Aporosa symplocoides", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo2, lwd = 2)
lines(needham_w_transform(K = S_parms["Elaeocarpus mastersii","K"], r = S_parms["Elaeocarpus mastersii", "r1"], p = S_parms["Elaeocarpus mastersii", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo2, lwd = 2)
mtext(side = 3, line = -1.3, text = " f)", adj = 0)
legend('bottomright', bty = "n", legend = c("High PC2: LMU, SCE", "Low PC2: ASY, EMA"), col = c(col.hi2, col.lo2), lwd = 2)

# g) Recruitment (Demog PC3)
b.bar <- traits_sp$rec[match(c("AFR", "GAX", "ACL", "MGI"),traits_sp$Species)]^2
barplot(b.bar, beside = T, col = c(col.hi3, col.hi3, col.lo3, col.lo3),
        ylab = expression(atop("BNR", "(recruits" ~ year^-1 ~ m^-2 ~ ")")), 
        xlab = "Species", ylim = c(0,25),
        names.arg = c("AFR", "GAX", "ACL", "MGI"), cex.lab = 1.5, cex.names = 1.2)
box()
mtext(side = 3, line = -1.3, text = " g)", adj = 0)

# h) Growth (Demog PC3)
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.8), xlab="DBH (cm)", ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Gynotroches axillaris","a"], b = AG_parms["Gynotroches axillaris", "b"], c = AG_parms["Gynotroches axillaris", "c"], dbh = newdbh) ~
        newdbh, col = col.hi3, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Aporosa frutescens","a"], b = AG_parms["Aporosa frutescens", "b"], c = AG_parms["Aporosa frutescens", "c"], dbh = newdbh) ~
        newdbh, col = col.hi3, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Archidendron clypearia","a"], b = AG_parms["Archidendron clypearia", "b"], c = AG_parms["Archidendron clypearia", "c"], dbh = newdbh) ~
        newdbh, col = col.lo3, lwd = 2)
lines(zeide_w_transform(a = AG_parms["Macaranga gigantea","a"], b = AG_parms["Macaranga gigantea", "b"], c = AG_parms["Macaranga gigantea", "c"], dbh = newdbh) ~ 
        newdbh, col = col.lo3, lwd = 2)
mtext(side = 3, line = -1.3, text = " h)", adj = 0)

# i) Survival (Demog PC3)
newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Macaranga gigantea","K"], r = S_parms["Macaranga gigantea", "r1"], p = S_parms["Macaranga gigantea", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo3, lwd = 2)
lines(needham_w_transform(K = S_parms["Archidendron clypearia","K"], r = S_parms["Archidendron clypearia", "r1"], p = S_parms["Archidendron clypearia", "p1"], dbh = newdbh) ~
        newdbh, col = col.lo3, lwd = 2)
lines(needham_w_transform(K = S_parms["Aporosa frutescens","K"], r = S_parms["Aporosa frutescens", "r1"], p = S_parms["Aporosa frutescens", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi3, lwd = 2)
lines(needham_w_transform(K = S_parms["Gynotroches axillaris","K"], r = S_parms["Gynotroches axillaris", "r1"], p = S_parms["Gynotroches axillaris", "p1"], dbh = newdbh) ~
        newdbh, col = col.hi3, lwd = 2)
mtext(side = 3, line = -1.3, text = " i)", adj = 0)
legend('bottomright', bty = "n", legend = c("High PC3: AFR, GAX", "Low PC3: ACL, MGI"), col = c(col.hi3, col.lo3), lwd = 2)

dev.off()

