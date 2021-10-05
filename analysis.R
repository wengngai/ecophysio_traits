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
pairwise <- data.frame(t(combn(ncol(traits_datonly),2)))
names(pairwise) <- c("y", "x")

for(i in 1:nrow(pairwise)){
    y <- traits_datonly[, pairwise[i,"y"] ]
    x <- traits_datonly[, pairwise[i,"x"] ]
    mod <- gls(y ~ x, correlation=corBrownian(phy = tree), method="ML")
    pairwise$y.name[i] <- names(traits_datonly)[ pairwise[i,"y"] ]
    pairwise$x.name[i] <- names(traits_datonly)[ pairwise[i,"x"] ]
    pairwise$coef[i] <- coef(mod)["x"]
    pairwise$p.value[i] <- summary(mod)$tTable["x", "p-value"]
    pairwise$AIC[i] <- summary(mod)$AIC
}
pairwise[which(pairwise$p.value < 0.05),]

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

newdbh <- seq(1,50,len=100)

# demog.PC1
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,1), xlab="DBH (log-transformed)", ylab="AGR", type="n")
lines(zeide_w_transform(a = AG_parms["Alstonia angustifolia","a"], b = AG_parms["Alstonia angustifolia", "b"], c = AG_parms["Alstonia angustifolia", "c"], dbh = newdbh) ~ 
          newdbh, col = "green")
lines(zeide_w_transform(a = AG_parms["Prunus polystachya","a"], b = AG_parms["Prunus polystachya", "b"], c = AG_parms["Prunus polystachya", "c"], dbh = newdbh) ~
          newdbh, col = "red")
lines(zeide_w_transform(a = AG_parms["Campnosperma squamatum","a"], b = AG_parms["Campnosperma squamatum", "b"], c = AG_parms["Campnosperma squamatum", "c"], dbh = newdbh) ~
          newdbh, col = "brown")
lines(zeide_w_transform(a = AG_parms["Rhodamnia cinerea","a"], b = AG_parms["Rhodamnia cinerea", "b"], c = AG_parms["Rhodamnia cinerea", "c"], dbh = newdbh) ~
          newdbh, col = "blue")
lines(zeide_w_transform(a = AG_parms["Macaranga bancana","a"], b = AG_parms["Macaranga bancana", "b"], c = AG_parms["Macaranga bancana", "c"], dbh = newdbh) ~
          newdbh, col = "grey")

# demog.PC2
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,1), xlab="DBH (log-transformed)", ylab="AGR", type="n")
lines(zeide_w_transform(a = AG_parms["Archidendron clypearia","a"], b = AG_parms["Archidendron clypearia", "b"], c = AG_parms["Archidendron clypearia", "c"], dbh = newdbh) ~
          newdbh, col = "brown")
lines(zeide_w_transform(a = AG_parms["Pternandra echinata","a"], b = AG_parms["Pternandra echinata", "b"], c = AG_parms["Pternandra echinata", "c"], dbh = newdbh) ~
          newdbh, col = "red")
lines(zeide_w_transform(a = AG_parms["Garcinia parvifolia","a"], b = AG_parms["Garcinia parvifolia", "b"], c = AG_parms["Garcinia parvifolia", "c"], dbh = newdbh) ~
          newdbh, col = "turquoise")
lines(zeide_w_transform(a = AG_parms["Lophopetalum multinervium","a"], b = AG_parms["Lophopetalum multinervium", "b"], c = AG_parms["Lophopetalum multinervium", "c"], dbh = newdbh) ~
          newdbh, col = "blue")

needham_w_transform <- function(K, r, p, dbh) K / (1 + exp(-r * (dbh - p) ))
newdbh <- seq(0.001, 1, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival probability", type="n")
lines(needham_w_transform(K = S_parms["Archidendron clypearia","K"], r = S_parms["Archidendron clypearia", "r1"], p = S_parms["Archidendron clypearia", "p1"], dbh = newdbh) ~
          newdbh, col = "brown")
lines(needham_w_transform(K = S_parms["Prunus polystachya","K"], r = S_parms["Prunus polystachya", "r1"], p = S_parms["Prunus polystachya", "p1"], dbh = newdbh) ~
          newdbh, col = "red")
lines(needham_w_transform(K = S_parms["Gynotroches axillaris","K"], r = S_parms["Gynotroches axillaris", "r1"], p = S_parms["Gynotroches axillaris", "p1"], dbh = newdbh) ~
          newdbh, col = "blue")
lines(needham_w_transform(K = S_parms["Aporosa symplocoides","K"], r = S_parms["Aporosa symplocoides", "r1"], p = S_parms["Aporosa symplocoides", "p1"], dbh = newdbh) ~
          newdbh, col = "turquoise")



#########################
# JSDMs with Phylo Dist #
#########################

traits_datonly <- apply(traits_datonly, 2, scale)

### HMSC (can do model selection) ###
library(Hmsc)



### Boral (cannot do model selection) ###

library(boral)

pd_mat <- cophenetic(tree)

# compare lv correlation structures

null <- boral(y = traits_datonly, X = cbind(demog.PC1, demog.PC2, demog.PC3, traits_sp$SSI), 
      family = "normal", calc.ics = T, 
      mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
      lv.control = list(num.lv = 2, type = "independent", distmat = NULL)
)

phylo.exp <- boral(y = traits_datonly, X = cbind(demog.PC1, demog.PC2, demog.PC3, traits_sp$SSI), 
                   family = "normal", calc.ics = T,
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "exponential", distmat = pd_mat)
)

phylo.sqexp <- boral(y = traits_datonly, X = cbind(demog.PC1, demog.PC2, demog.PC3, traits_sp$SSI), 
                   family = "normal", 
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "squared.exponential", distmat = pd_mat)
)

phylo.powexp <- boral(y = traits_datonly, X = cbind(demog.PC1, demog.PC2, demog.PC3, traits_sp$SSI), 
                   family = "normal", 
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "powered.exponential", distmat = pd_mat)
)

phylo.spher <- boral(y = traits_datonly, X = cbind(demog.PC1, demog.PC2, demog.PC3, traits_sp$SSI), 
                   family = "normal", 
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "spherical", distmat = pd_mat)
)


offset <- c(-0.3, -0.15, 0, 0.15, 0.3)
library(viridis)
cols <- viridis(5)
colt <- adjustcolor(cols, alpha.f = 0.25)
nt <- ncol(traits_datonly)

pdf("D:/Dropbox/Functional Traits Project/Figures/HMSC coef plot Sep21.pdf", height = 9, width = 16)
par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(1,7,1,1))
for(i in 1:4){
  plot(x = 1:nt, y = 1:nt, type = "n",
       xlim = c(-4,4), xlab = "", yaxt = "n", ylab = "")
  abline(v = 0, lty = 2)
  for(j in 1:5){
    if(j==1) mod <- null
    if(j==2) mod <- phylo.exp
    if(j==3) mod <- phylo.sqexp
    if(j==4) mod <- phylo.powexp
    if(j==5) mod <- phylo.spher
    colsig <- rep(colt[j], nt)
    sig.index <- which(mod$hpdintervals$X.coefs[,i,"lower"] > 0 | mod$hpdintervals$X.coefs[,i,"upper"] < 0)
    colsig[sig.index] <- cols[j]
    ##colsig <- which(sig==1)
    points(c(1:nt) + offset[j] ~ mod$X.coefs.median[,i], pch = 16, cex = 2, col = colsig)
    for(k in 1:nt){
      arrows(mod$hpdintervals$X.coefs[k,i,"lower"], k + offset[j], mod$hpdintervals$X.coefs[k,i,"upper"], k + offset[j], 
             length = 0, col = colsig[k])
    }
  }
  if(i==1){
    mtext(side=3, text = "a) Demographic PC1", cex = 1.5, adj = 0)
    axis(side=2, at = 1:nt, labels = colnames(traits_datonly), las = 1)
  }
  if(i==2) mtext(side=3, text = "b) Demographic PC2", cex = 1.5, adj = 0)
  if(i==3){
    mtext(side=3, text = "c) Demographic PC3", cex = 1.5, adj = 0)
    axis(side=2, at = 1:nt, labels = colnames(traits_datonly), las = 1)
  }
  if(i==4) mtext(side=3, text = "d) SSI", cex = 1.5, adj = 0)
}
legend('bottomright', bty = "n", pch = 16, col = rev(cols), pt.cex = 2, legend = rev(c(
  "No phylogenetic structure",
  "Exponential",
  "Squared exponential",
  "Powered exponential",
  "Spherical"))
)
mtext(side = 2, outer = T, "Trait", line = 4, cex = 1.5)
mtext(side = 1, outer = T, "Effect size", line = -1, cex = 1.5)







# Clustering traits?
traits.clust <- hclust(dist(traits_datonly))
plot(traits.clust, labels = traits_sp$Species)




