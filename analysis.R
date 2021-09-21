library(vegan)
traits_sp <- read.csv("combined traits_sp level_Sep21.csv", row.names = 1, header = T)
summary(traits_sp)


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
#tree$tip.label <- abbrev.tip(tree$tip.label)
#tree$tip.label <- spp_list$Species[match(abbrev.tip(tree$tip.label), spp_list$Sp)]
tree$tip.label <- spp_list$long_name[match(abbrev.tip(tree$tip.label), spp_list$Sp)]

plot(tree, no.margin=TRUE)


### Demographic params PCA ###
PCA.demog <- rda(traits_sp[26:31], scale=T)
summary(PCA.demog)$cont
biplot(PCA.demog, choices=c(1,2))
biplot(PCA.demog, choices=c(2,3))

demog.PC1 <- PCA.demog$CA$u[,1]
demog.PC2 <- PCA.demog$CA$u[,2]
demog.PC3 <- PCA.demog$CA$u[,3]

### Traits PCA ###
PCA.traits <- rda(traits_sp[2:25], scale = T)
summary(PCA.traits)$cont
rownames(PCA.traits$CA$u) <- traits_sp$Species
biplot(PCA.traits, choices = c(1,2))
biplot(PCA.traits, choices = c(2,3))

# RDA of traits against SSI and demographic params
rda.ssi <- rda(traits_sp[2:25] ~ traits_sp$SSI + demog.PC1 + demog.PC2 + demog.PC3, scale = T)
anova(rda.ssi, by="margin")

### Testing all possible trait-trait correlations ###
library(nlme)
library(ape)

pairwise <- data.frame(t(combn(ncol(traits_sp)-1,2)))
names(pairwise) <- c("y", "x")
traits_datonly <- traits_sp[,-1]

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


# GLS of demog params against trait PCs

tPC1 <- PCA.traits$CA$u[,1]
tPC2 <- PCA.traits$CA$u[,2]
tPC3 <- PCA.traits$CA$u[,3]
tPC4 <- PCA.traits$CA$u[,4]
tPC5 <- PCA.traits$CA$u[,5]
tPC6 <- PCA.traits$CA$u[,6]

summary(gls(demog.PC1 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6, correlation=corBrownian(phy = tree), method="ML"))
summary(gls(demog.PC2 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6, correlation=corBrownian(phy = tree), method="ML"))
summary(gls(demog.PC3 ~ tPC1 + tPC2 + tPC3 + tPC4 + tPC5 + tPC6, correlation=corBrownian(phy = tree), method="ML"))

