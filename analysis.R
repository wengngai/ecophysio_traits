### PCA ###
PCA.1 <- prcomp(traits_sp[2:length(traits_sp)], scale=T)
biplot(PCA.1)
summary(PCA.1)

library(vegan)
rda.ssi <- rda(traits_sp[2:22] ~ traits_sp$SSI)
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
