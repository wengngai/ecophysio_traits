library(vegan)
traits_sp <- read.csv("combined traits_sp level_Sep21.csv", row.names = 1, header = T)
summary(traits_sp)

traits_sp[c("Species", "SSI")]

### Phylogenetic Tree ###

library(brranching)

#Input taxonomic names

phylolist <- c(
    "Alstonia angustifolia",
    "Aporosa frutescens",
    "Aporosa symplocoides",
    "Pithecellobium clypearia",
    "Baccaurea bracteata",
    "Bhesa paniculata",
    "Campnosperma squamatum",
    "Elaeocarpus mastersii",
    "Garcinia parvifolia",
    "Gironniera nervosa",
    "Gynotroches axillaris",
    "Horsfieldia crassifolia",
    "Knema malayana",
    "Lophopetalum multinervium",
    "Macaranga bancana",
    "Macaranga gigantea",
    "Macaranga recurvata",
    "Mussaendopsis beccariana",
    "Pellacalyx axillaris",
    "Pometia pinnata",
    "Prunus polystachya",
    "Pternandra coerulescens",
    "Pternandra echinata",
    "Rhodamnia cinerea",
    "Strombosia ceylanica",
    "Syzygium pachyphyllum",
    "Timonius flavescens",
    "Timonius wallichianus",
    "Xanthophyllum flavescens"
)

tree = phylomatic(taxa=phylolist, storedtree = "zanne2014")

# rename tips
spp_list <- read.csv("./raw_data/species list.csv")
spp_list$long_name <- paste0(spp_list$Species, " (", spp_list$Family, ")")

library(stringr)
abbrev.tip <- function(x){
    names <- str_split(x, "_")
    abbrev.phylo <- toupper(paste0(substr(sapply(names, "[[", 1), 1, 1), substr(sapply(names, "[[", 2), 1, 2)))
    abbrev.phylo[which(abbrev.phylo=="PCL")] <- "ACL"
    return(abbrev.phylo)
}

# sort traits_sp into order of phylo tree first
traits_sp <- traits_sp[match(abbrev.tip(tree$tip.label), traits_sp$Species),]

# EITHER OR:
tree$tip.label <- abbrev.tip(tree$tip.label)
#tree$tip.label <- spp_list$Species[match(abbrev.tip(tree$tip.label), spp_list$Sp)]
#tree$tip.label <- spp_list$long_name[match(abbrev.tip(tree$tip.label), spp_list$Sp)]

#pdf("./outputs/Phylo tree.pdf", width=8, height=8)
plot(tree, no.margin=TRUE, 
     tip.color = ifelse(traits_sp$SSI > 0.66, "steelblue", ifelse(traits_sp$SSI < 0.33, "forestgreen", "grey40")))
dev.off()



### Demographic params PCA ###
# include Hmax as demographic param or not? it doesn't really add much
PCA.demog <- rda(traits_sp[26:31], scale=T)
summary(PCA.demog)$cont
rownames(PCA.demog$CA$u) <- traits_sp$Species

#pdf("./outputs/demog PCA1-2.pdf", width=5, height=5)
biplot(PCA.demog, choices=c(1,2), cex=3, cex.lab=1.5) 
dev.off()
# PC1 = growth rate (higher PC1 = faster growth rates)
# PC2 = attenuation of growth rate in larger trees (higher PC2 = less attenuation)
#pdf("./outputs/demog PCA1-3.pdf", width=5, height=5)
biplot(PCA.demog, choices=c(1,3))
dev.off()
# PC3 = under- (neg) vs over- (pos) canopy species (growth-survival tradeoff axis)
# Even if Hmax is included, it does NOT load on to PC3, even though many subcanopy species have negative PC3

demog.PC1 <- PCA.demog$CA$u[,1]
demog.PC2 <- PCA.demog$CA$u[,2]
demog.PC3 <- PCA.demog$CA$u[,3]


### Traits PCA ###
PCA.traits <- rda(traits_sp[2:25], scale = T)
summary(PCA.traits)$cont
rownames(PCA.traits$CA$u) <- traits_sp$Species

#pdf("./outputs/traits PCA1-2.pdf", width=5, height=5)
biplot(PCA.traits, choices = c(1,2))
dev.off()
# PC1: acquisitive-conservative spectrum
# PC2: hydraulic safety versus 
biplot(PCA.traits, choices = c(1, 4))

# RDA of traits against SSI and demographic params
rda.ssi <- rda(traits_sp[2:25] ~ traits_sp$SSI + demog.PC1 + demog.PC2 + demog.PC3, scale = T)
anova(rda.ssi, by="margin")

### Exploring all possible trait-trait correlations ###
library(nlme)
library(ape)

traits_datonly <- traits_sp[,2:25]
traits_datonly_scaled <- apply(traits_datonly, 2, scale)
pairwise <- data.frame(t(combn(ncol(traits_datonly_scaled),2)))
names(pairwise) <- c("y", "x")

for(i in 1:nrow(pairwise)){
    y <- traits_datonly_scaled[, pairwise[i,"y"] ]
    x <- traits_datonly_scaled[, pairwise[i,"x"] ]
    mod <- gls(y ~ x, correlation=corBrownian(phy = tree), method="ML")
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
  "WD", "Hmax"
)
groupnames <- c(rep("leaf", 13), rep("wood", 11))
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

# Fit GLS models
library(MuMIn)

# demog.PC1
summary(max.dPC1 <- gls(demog.PC1 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
            correlation=corBrownian(phy = tree), method="ML"))
(dr.dPC1 <- dredge(max.dPC1))
summary(best.dPC1 <- get.models(dr.dPC1, subset=1)[[1]])
#pdf("./outputs/traits PCA2-4.pdf", width=5, height=5)
biplot(PCA.traits, choices=c(2,4), display = "species") 
dev.off()

# demog.PC2
summary(max.dPC2 <- gls(demog.PC2 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
            correlation=corBrownian(phy = tree), method="ML"))
(dr.dPC2 <- dredge(max.dPC2))
summary(best.dPC2 <- get.models(dr.dPC2, subset=1)[[1]])

# demog.PC3
summary(max.dPC3 <- gls(demog.PC3 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
            correlation=corBrownian(phy = tree), method="ML"))
(dr.dPC3 <- dredge(max.dPC3))
summary(best.dPC3 <- get.models(dr.dPC3, subset=1)[[1]])

# SSI
summary(max.ssi <- gls(model = SSI ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6 + tPC7 + tPC8, 
            data = traits_sp, correlation=corBrownian(phy = tree), method="ML"))
(dr.ssi <- dredge(max.ssi))
summary(best.ssi <- get.models(dr.ssi, subset=1)[[1]])

#pdf("./outputs/traits PCA 4-7.pdf", width=5, height=5)
biplot(PCA.traits, choices=c(7,4), display = "species") 
dev.off()



# Plot the significant correlations
plot(demog.PC1 ~ tPC2, type="n"); text(demog.PC1 ~ tPC2, labels=traits_sp$Species)
plot(demog.PC1 ~ tPC4, type="n"); text(demog.PC1 ~ tPC4, labels=traits_sp$Species)
plot(demog.PC2 ~ tPC1, type="n"); text(demog.PC2 ~ tPC1, labels=traits_sp$Species)
plot(demog.PC3 ~ tPC1, type="n"); text(demog.PC3 ~ tPC1, labels=traits_sp$Species)
plot(demog.PC3 ~ tPC2, type="n"); text(demog.PC3 ~ tPC2, labels=traits_sp$Species)
plot(demog.PC3 ~ tPC3, type="n"); text(demog.PC3 ~ tPC3, labels=traits_sp$Species)
plot(demog.PC3 ~ tPC4, type="n"); text(demog.PC3 ~ tPC4, labels=traits_sp$Species)
plot(traits_sp$SSI ~ tPC1, type="n"); text(traits_sp$SSI ~ tPC1, labels=traits_sp$Species)
plot(traits_sp$SSI ~ tPC4, type="n"); text(traits_sp$SSI ~ tPC4, labels=traits_sp$Species)
plot(traits_sp$SSI ~ tPC7, type="n"); text(traits_sp$SSI ~ tPC7, labels=traits_sp$Species)

biplot(PCA.traits, choices = c(4,7))


# are demog.PCs and SSI correlated?
summary(ssi.mod <- gls(model = SSI ~ demog.PC1 + demog.PC2 + demog.PC3,
                       data = traits_sp, correlation=corBrownian(phy = tree), method="ML"))
dr.ssi.demo <- dredge(ssi.mod)
summary(get.models(dr.ssi.demo, subset=1)[[1]])


#################
# VISUALIZATION #
#################

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
needham_w_transform <- function(K, r, p, dbh) K / (1 + exp(-r * (dbh - p) ))

### PCA plot
#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Demog PCA.pdf", height=6, width=8)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Demog PCA.jpg", height=6, width=8, units="in", res=300)
layout(matrix(c(1,1,2,1,1,3,4,5,6), ncol = 3))
par(mar = c(5,5.5,2,2))
plot(demog.PC2 ~ demog.PC1, type = "n", ylim = c(-0.45,0.45), xlim = c(-0.6, 0.4),
     ylab = "PC2 (20.5% vriance explained)", xlab = "PC1 (43.8% variance explained)", cex.lab = 1.5)
for(i in 1:6){
  if(i<=3) arrowcol <- "#403D58" else arrowcol <- "#8E8963" 
  arrows(0, 0, PCA.demog$CA$v[i,1]/2, PCA.demog$CA$v[i,2]/2, lwd = 2, col = arrowcol)
  text(x = PCA.demog$CA$v[i,1]/1.8, y = PCA.demog$CA$v[i,2]/1.8, lwd = 2, col = arrowcol, 
       labels = rownames(PCA.demog$CA$v)[i], cex = 1.2 )
}
text(demog.PC2 ~ demog.PC1, labels = traits_sp$Species,
     col = ifelse(traits_sp$Species %in% c("LMU","GPA"), "#F7B39F", 
                  ifelse(traits_sp$Species %in% c("PEC","EMA"), "#66D7D1",
                         ifelse(traits_sp$Species %in% c("PPO", "GNE"), "#FC7753",
                                ifelse(traits_sp$Species %in% c("RCI", "AAN"), "#538A95", "grey")))),
     cex = ifelse(traits_sp$Species %in% c("LMU", "GPA", "PEC", "EMA", "PPO", "GNE", "RCI", "AAN"), 1, 0.7)
     )
legend("bottomleft", bty = "n", lwd = 2, col=c("#403D58","#8E8963"),
       legend = c("Growth model parameters: a, b, c (panels b & d)", "Survival model parameters: K, r1, p1 (panels c & e)"))
mtext(side = 3, line = -1.2, text = " a)", adj = 0)

## Demog PC1
# Growth
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.7), xlab="DBH (cm)", ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Alstonia angustifolia","a"], b = AG_parms["Alstonia angustifolia", "b"], c = AG_parms["Alstonia angustifolia", "c"], dbh = newdbh) ~ 
        newdbh, col = "#538A95", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Prunus polystachya","a"], b = AG_parms["Prunus polystachya", "b"], c = AG_parms["Prunus polystachya", "c"], dbh = newdbh) ~
        newdbh, col = "#FC7753", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Gironniera nervosa","a"], b = AG_parms["Gironniera nervosa", "b"], c = AG_parms["Gironniera nervosa", "c"], dbh = newdbh) ~
        newdbh, col = "#FC7753", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Rhodamnia cinerea","a"], b = AG_parms["Rhodamnia cinerea", "b"], c = AG_parms["Rhodamnia cinerea", "c"], dbh = newdbh) ~
        newdbh, col = "#538A95", lwd = 2)
#text(x = rep(20,4), y = c(0.65, 0.48, 0.26, 0.17), labels = c("PPO", "GNE", "AAN", "RCI"), col=c("#FC7753", "#FC7753", "#538A95", "#538A95"))
mtext(side = 3, line = -1.2, text = " b)", adj = 0)

# Survival
newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival probability", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Gironniera nervosa","K"], r = S_parms["Gironniera nervosa", "r1"], p = S_parms["Gironniera nervosa", "p1"], dbh = newdbh) ~
        newdbh, col = "#FC7753", lwd = 2)
lines(needham_w_transform(K = S_parms["Prunus polystachya","K"], r = S_parms["Prunus polystachya", "r1"], p = S_parms["Prunus polystachya", "p1"], dbh = newdbh) ~
        newdbh, col = "#FC7753", lwd = 2)
lines(needham_w_transform(K = S_parms["Alstonia angustifolia","K"], r = S_parms["Alstonia angustifolia", "r1"], p = S_parms["Alstonia angustifolia", "p1"], dbh = newdbh) ~
        newdbh, col = "#538A95", lwd = 2)
lines(needham_w_transform(K = S_parms["Rhodamnia cinerea","K"], r = S_parms["Rhodamnia cinerea", "r1"], p = S_parms["Rhodamnia cinerea", "p1"], dbh = newdbh) ~
        newdbh, col = "#538A95", lwd = 2)
legend('bottomright', bty = "n", legend = c("Low PC1: AAN, RCI", "High PC1: GNE, PPO"), col = c("#538A95", "#FC7753"), lwd = 2)
mtext(side = 3, line = -1.2, text = " c)", adj = 0)

## Demog PC2
# Growth
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.45), xlab="DBH (cm)", ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Lophopetalum multinervium","a"], b = AG_parms["Lophopetalum multinervium", "b"], c = AG_parms["Lophopetalum multinervium", "c"], dbh = newdbh) ~ 
        newdbh, col = "#F7B39F", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Elaeocarpus mastersii","a"], b = AG_parms["Elaeocarpus mastersii", "b"], c = AG_parms["Elaeocarpus mastersii", "c"], dbh = newdbh) ~
        newdbh, col = "#66D7D1", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Garcinia parvifolia","a"], b = AG_parms["Garcinia parvifolia", "b"], c = AG_parms["Garcinia parvifolia", "c"], dbh = newdbh) ~ 
        newdbh, col = "#F7B39F", lwd = 2)
lines(zeide_w_transform(a = AG_parms["Pternandra echinata","a"], b = AG_parms["Pternandra echinata", "b"], c = AG_parms["Pternandra echinata", "c"], dbh = newdbh) ~
        newdbh, col = "#66D7D1", lwd = 2)
mtext(side = 3, line = -1.2, text = " d)", adj = 0)

# Survival
newdbh <- seq(0.001, 0.5, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival probability", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Lophopetalum multinervium","K"], r = S_parms["Lophopetalum multinervium", "r1"], p = S_parms["Lophopetalum multinervium", "p1"], dbh = newdbh) ~
        newdbh, col = "#F7B39F", lwd = 2)
lines(needham_w_transform(K = S_parms["Elaeocarpus mastersii","K"], r = S_parms["Elaeocarpus mastersii", "r1"], p = S_parms["Elaeocarpus mastersii", "p1"], dbh = newdbh) ~
        newdbh, col = "#66D7D1", lwd = 2)
lines(needham_w_transform(K = S_parms["Garcinia parvifolia","K"], r = S_parms["Garcinia parvifolia", "r1"], p = S_parms["Garcinia parvifolia", "p1"], dbh = newdbh) ~
        newdbh, col = "#F7B39F", lwd = 2)
lines(needham_w_transform(K = S_parms["Pternandra echinata","K"], r = S_parms["Pternandra echinata", "r1"], p = S_parms["Pternandra echinata", "p1"], dbh = newdbh) ~
        newdbh, col = "#66D7D1", lwd = 2)
mtext(side = 3, line = -1.2, text = " e)", adj = 0)
legend('bottomright', bty = "n", legend = c("Low PC2: EMA, PEC", "High PC2: GPA, LMU"), col = c("#66D7D1", "#F7B39F"), lwd = 2)

dev.off()
