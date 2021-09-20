library(vegan)
traits_sp <- read.csv("combined traits_sp level_Sep21.csv", row.names = 1, header = T)
names(traits_sp)
summary(traits_sp)

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
biplot(PCA.traits, choices = c(1,2))

rda.ssi <- rda(traits_sp[2:25] ~ traits_sp$SSI + demog.PC1 + demog.PC2 + demog.PC3, scale = T)
anova(rda.ssi, by="margin")

### Phylogenetic Tree ###

library(brranching)

#Input taxonomic names

spplist <- data.frame(
    full = c(
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
    ),
    abbr = c(
        "AAN",
        "AFR",
        "ASY",
        "ACL",
        "BBR",
        "BPA",
        "CSQ",
        "EMA",
        "GPA",
        "GNE",
        "GAX",
        "HCR",
        "KMA",
        "LMU",
        "MBA",
        "MGI",
        "MRE",
        "MBE",
        "PAX",
        "PPI",
        "PPO",
        "PCO",
        "PEC",
        "RCI",
        "SCE",
        "SPA",
        "TFL",
        "TWA",
        "XFL"
    )
)

tree = phylomatic(taxa=spplist$full, storedtree = "zanne2014")

plot(tree, no.margin=TRUE)
