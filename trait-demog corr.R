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

### Traits ###

names(traits_sp)
traits_sp[,c(15,18:23)]

############################
# TRAIT-TRAIT CORRELATIONS #
############################

Tz <- traits_sp[2:25]
traits_datonly <- Tz
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
pairwise$p.value <- p.adjust(pairwise$p.value, method = "fdr", choose(ncol(Tz), 2))
pairwise[which(pairwise$p.value < 0.05),]

library(circlize)

ordered.traitnames <- c(
    "SLA", "LDMC", "CNR",
    "Th", "Th_LC", "Th_LE", "Th_SM", "Th_PM", "Th_UE", "Th_UC", "L_VD",
    "GCL", "SD", 
    "RV", "CV", "SV", "OV", "TV", "VGI", "VA", "DH", "T_VD",
    "WD", "Hmax"
)
groupnames <- c(rep("leaf", 13), rep("wood", 11))
grid.col <- ifelse(groupnames=="leaf", "#DBD56E", "#403D58")
names(groupnames) <- ordered.traitnames
names(grid.col) <- ordered.traitnames

pairwise$coef[which(pairwise$p.value > 0.05)] <- 0
links.col <- ifelse(pairwise$coef>0, "#66D7D1", ifelse(pairwise$coef==0, "#F2EFEA", "#FC7753"))
alpha.sig <- (0.05 - pairwise$p.value)/0.1
alpha.sig[which(alpha.sig < 0)] <- 0
for(i in 1:length(links.col)) links.col[i] <- adjustcolor(links.col[i], alpha.f = alpha.sig[i])

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait correlations chord FDR.pdf", height=6, width=10)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Trait correlations chord FDR.jpg", height=6, width=10, units="in", res=300)
chordDiagram(pairwise[,3:5], symmetric = T, 
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

lat <- ordered.traitnames[5:10]

lat.pw <- pairwise[which(pairwise$y.name %in% lat & pairwise$x.name %in% lat), 3:5]

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Leaf layer correlations.pdf", height=6, width=10)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Leaf layer correlations.jpg", height=6, width=10, units="in", res=300)
chordDiagram(lat.pw,
             symmetric = T, order = lat,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = mm_h(c(3, 3)),
             grid.col = "#DBD56E", transparency = 0.5, 
             col = ifelse(lat.pw$coef>0, "#66D7D1", ifelse(lat.pw$coef==0, "#F2EFEA", "#FC7753"))
)
legend(x= 1.1, y = -0.5, title = "Correlation", bty = "n",
       legend=c("Positive", "Negative"), fill = adjustcolor(c("#66D7D1", "#FC7753"), alpha.f = 0.5))
dev.off()
circos.clear()


pairwise[pairwise$p.value < 0.05, 3:5]


sel <- c("SLA", "LDMC", "CNR", "Hmax", "WD",
         "GCL", "SD", "Th_SM", "Th", "L_VD",
         "T_VD", "VGI", "VA"
)

########
# PGLS #
########

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

demog <- cbind(traits_sp[26:32], traits_sp$logitSSI)
Z <- traits_sp[sel]
output <- vector(mode = 'list', length = ncol(demog))

for(i in 1:ncol(demog)){
    output[[i]] <- data.frame(coef = rep(NA, ncol(Z)), upp = NA, low = NA, y = ncol(Z):1)
    rownames(output[[i]]) <- names(Z)
    for(j in 1:ncol(Z)){
        Y <- scale(demog[,i])
        X <- scale(Z[j])
        mod <- gls(Y ~ X, correlation=corBrownian(phy = tree$scenario.3), method="ML")
        output[[i]]$coef[j] <- summary(mod)$tTable[2,1]
        output[[i]]$upp[j] <- summary(mod)$tTable[2,1] + 1.96*summary(mod)$tTable[2,2]
        output[[i]]$low[j] <- summary(mod)$tTable[2,1] - 1.96*summary(mod)$tTable[2,2]
    }
}

names(demog)
panellabels <- c("a) Growth: a",
                "b) Growth: b",
                "c) Growth: c",
                "d) Survival: K",
                "e) Survival: r",
                "f) Survival: p",
                "g) Sapling recruitment: BNR",
                "h) Swamp specialization index")
colmat <- matrix(c("darkgreen", "greenyellow", "darkgreen", "greenyellow", "darkgreen", "greenyellow",
            "darkred", "lightpink", "darkred", "lightpink", "darkred", "lightpink",
            "gold4", "khaki", "darkblue", "lightblue"), nrow = 8, ncol = 2, byrow = T)

#pdf("./outputs/trait-demog correlations.pdf", height=10, width=6)
#jpeg("./outputs/trait-demog correlations.jpg", height=10, width=6, units = "in", res = 300)
par(mfrow = c(4,2), mar = c(3,2,2,2), oma = c(2,3,1,6))
for(i in 1:ncol(demog)){
    cols <- ifelse(output[[i]]$upp < 0 | output[[i]]$low > 0, colmat[i,1], colmat[i,2])
    if(i %in% c(1,3,5,7)) par(mar = c(3,3,2,1)) else par(mar = c(3,1,2,3))
        plot(y ~ coef, pch = 16, xlim = c(min(low),max(upp)), data = output[[i]], col = cols, 
         yaxt = "n", xlab = "", cex = 2)
    if(i %in% c(1,3,5,7)) axis(side = 2, at = output[[i]]$y, labels = rownames(output[[i]]), las = 1) else{
        axis(side = 4, at = c(2,5,7.5,9.5,12), las = 1, tick = F, 
             labels = c("Twig anatomical", "Leaf anatomical", "Stomatal", "Wood", "LES"))
        axis(side = 2, at = output[[i]]$y, labels = NA, las = 1)
        axis(side = 4, at = c(3.5, 6.5, 8.5, 10.5), labels = NA)
    }
    for(j in 1:ncol(Z)){
        arrows(output[[i]]$low[j], output[[i]]$y[j], output[[i]]$upp[j], output[[i]]$y[j],
               length = 0, col = cols[j], lwd = 2)
    }
    abline(v = 0, lty = 2)
    abline(h = c(3.5, 6.5, 8.5, 10.5), lty = 2)
    mtext(side = 3, adj = 0, text = panellabels[i])
}
mtext(side = 1, outer = T, text = "Standardized effect size")
dev.off()





######################
# PREDICTION FIGURES #
######################
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


# assign colors
col.hi1 <- "#66D7D1"
col.lo1 <- "#B1A792"
col.hi2 <- "#538A95"
col.lo2 <- "#FC7753"
col.hi3 <- "#403D58"
col.lo3 <- "#F7B39F"

#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\Zeide and Needham correlations.pdf", height=9, width=4.5)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\Zeide and Needham correlations.jpg", height=9, width=4.5, units = "in", res = 600)
par(mfrow = c(2,1), mar = c(4.5,5.5,1,1))
newdbh <- seq(1,50,len=100)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0,0.7), xlab="DBH (cm)", 
     ylab=expression(paste("AGR (cm ", year^-1, ")")), type="n", cex.lab = 1.5)
lines(zeide_w_transform(a = AG_parms["Rhodamnia cinerea","a"], b = AG_parms["Rhodamnia cinerea", "b"], c = AG_parms["Rhodamnia cinerea", "c"], dbh = newdbh) ~
          newdbh, col = col.lo3, lwd = 3)
lines(zeide_w_transform(a = AG_parms["Mussaendopsis beccariana","a"], b = AG_parms["Mussaendopsis beccariana", "b"], c = AG_parms["Mussaendopsis beccariana", "c"], dbh = newdbh) ~
          newdbh, col = col.lo3, lwd = 3, lty = 2)
lines(zeide_w_transform(a = AG_parms["Macaranga gigantea","a"], b = AG_parms["Macaranga gigantea", "b"], c = AG_parms["Macaranga gigantea", "c"], dbh = newdbh) ~
          newdbh, col = col.hi3, lwd = 3)
lines(zeide_w_transform(a = AG_parms["Aporosa symplocoides","a"], b = AG_parms["Aporosa symplocoides", "b"], c = AG_parms["Aporosa symplocoides", "c"], dbh = newdbh) ~
          newdbh, col = col.hi3, lwd = 3, lty = 2)
legend('topright', bty = "n",
       legend = c("High VA, high VGI: MGI", "Low VA, low VGI: RCI", "Low SD: ASY", "High SD: MBE"),
       lwd = 3, col = c(col.hi3, col.lo3), lty = c(1,1,2,2))
mtext(side = 3, adj = 0, line = -1.5, text = " a)", cex = 1.5)

newdbh <- seq(0.001, 1, len=200)
plot(rep(0.5, length(newdbh)) ~ newdbh, ylim=c(0.4,1), xlab="DBH (cm)", ylab="Survival", type="n", cex.lab = 1.5)
lines(needham_w_transform(K = S_parms["Aporosa symplocoides","K"], r = S_parms["Aporosa symplocoides", "r1"], p = S_parms["Aporosa symplocoides", "p1"], dbh = newdbh) ~
          newdbh, col = col.hi2, lwd = 3)
lines(needham_w_transform(K = S_parms["Campnosperma squamatum","K"], r = S_parms["Campnosperma squamatum", "r1"], p = S_parms["Campnosperma squamatum", "p1"], dbh = newdbh) ~
          newdbh, col = col.hi2, lwd = 3, lty = 2)
lines(needham_w_transform(K = S_parms["Archidendron clypearia","K"], r = S_parms["Archidendron clypearia", "r1"], p = S_parms["Archidendron clypearia", "p1"], dbh = newdbh) ~
          newdbh, col = col.lo2, lwd = 3, lty = 2)
lines(needham_w_transform(K = S_parms["Mussaendopsis beccariana","K"], r = S_parms["Mussaendopsis beccariana", "r1"], p = S_parms["Mussaendopsis beccariana", "p1"], dbh = newdbh) ~
          newdbh, col = col.lo2, lwd = 3)
legend('bottomright', bty = "n", legend = c("High SD: MBE", "Low SD: ASY", "High Th_SM: CSQ", "Low Th_SM: MBE"), col = c(col.hi2, col.lo2), lwd = 3, lty = c(1,1,2,2))
mtext(side = 3, line = -1.5, text = " d)", adj = 0, cex = 1.5)
# note: panels b and c are photos of twig cross sections
# panels e and f, of stomata
dev.off()





plot(rec ~ VGI, data = traits_sp, type = "n"); text(rec ~ VGI, data = traits_sp, labels = Species)
plot(VA ~ VGI, data = traits_sp, type = "n"); text(VA ~ VGI, data = traits_sp, labels = Species)
plot(SD ~ VGI, data = traits_sp, type = "n"); text(SD ~ VGI, data = traits_sp, labels = Species)
pairs.cor(traits_sp[c("SD", "VGI", "VA")])
plot(a ~ VGI, data = traits_sp, type = "n"); text(a ~ VGI, data = traits_sp, labels = Species)
plot(WD ~ Hmax, data = traits_sp, type = "n"); text(WD ~ Hmax, data = traits_sp, labels = Species)
plot(WD ~ K, data = traits_sp, type = "n"); text(WD ~ K, data = traits_sp, labels = Species)
plot(Hmax ~ K, data = traits_sp, type = "n"); text(Hmax ~ K, data = traits_sp, labels = Species)
plot(SD ~ K, data = traits_sp, type = "n"); text(SD ~ K, data = traits_sp, labels = Species)

plot(SSI ~ Th, data = traits_sp, type = "n"); text(SSI ~ Th, data = traits_sp, labels = Species)
plot(SSI ~ Th_SM, data = traits_sp, type = "n"); text(SSI ~ Th_SM, data = traits_sp, labels = Species)
plot(K ~ Th_SM, data = traits_sp, type = "n"); text(K ~ Th_SM, data = traits_sp, labels = Species)

plot(SD ~ rec, data = traits_sp, type = "n"); text(SD ~ rec, data = traits_sp, labels = Species)
plot(L_VD ~ rec, data = traits_sp, type = "n"); text(L_VD ~ rec, data = traits_sp, labels = Species)


hist(traits_sp$VGI)




# SSI #

library(rr2)
summary(SM.ssi <- gls(model = logitSSI ~ Th_SM, data = traits_sp, correlation=corBrownian(phy = tree$scenario.3), method="ML"))
# create null model to estimate R2
m0.ssi <- gls(logitSSI ~ 1, correlation=corBrownian(phy = tree$scenario.3), data = traits_sp, method="ML")
R2.lik(SM.ssi, m0.ssi)
# predictions
newSM <- seq(-0.8, 0.8, 0.01)
SMpred <- predict(SM.ssi, newdata = data.frame(Th_SM = newSM), se.fit = T)
inv.logit <- function(y) exp(y) / (1 + exp(y))

# Sapling recruitment #

summary(LVD.BSR <- gls(model = rec ~ L_VD, data = traits_sp, correlation=corBrownian(phy = tree$scenario.3), method="ML"))
# create null model to estimate R2
m0.BSR <- gls(rec ~ 1, correlation=corBrownian(phy = tree$scenario.3), data = traits_sp, method="ML")
R2.lik(LVD.BSR, m0.BSR)
# predictions
newLVD <- seq(0.005, 0.012, 0.0001)
LVDpred <- predict(LVD.BSR, newdata = data.frame(L_VD = newLVD), se.fit = T)


#pdf("D:\\Dropbox\\Functional Traits Project\\Figures\\SSI and BSR correlations.pdf", height=9, width=4.5)
#jpeg("D:\\Dropbox\\Functional Traits Project\\Figures\\SSI and BSR correlations.jpg", height=9, width=4.5, units = "in", res = 600)
par(mfrow = c(2,1), mar = c(4.5,5.5,1,1), mgp = c(3.5,1,0))
plot(SSI ~ Th_SM, type="n", data = traits_sp, cex.lab = 1.5,
     ylab = "Swamp specialization index", 
     xlab = "Spongy mesophyll\nrelative thickness (Th_SM)")
polygon(with(SMpred, inv.logit(c(fit + se.fit, rev(fit - se.fit)))) ~ c(newSM,rev(newSM)), border = F, col = "#66D7D16e")
lines(inv.logit(SMpred$fit) ~ newSM, lwd = 4, col = "white")
text(SSI ~ Th_SM, labels = Species, data = traits_sp,
     cex = ifelse(traits_sp$Species %in% c("CSQ", "ACL"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("CSQ", "ACL"), "black", "grey60"))
mtext(side = 3, adj = 0, line = -1.5, text = " a)", cex = 1.5)

plot(rec ~ L_VD, type="n", data = traits_sp, xlim = c(0.005,0.012), cex.lab = 1.5,
     ylab = expression(paste("BSR (stems ", year^-1, " ", plot^-1, " ", m^-2, ")")),
     xlab = expression(paste("Leaf vein density (", mu, "m ", mu, "m" ^-2, ")")))
polygon(with(LVDpred, c(fit + se.fit, rev(fit - se.fit))) ~ c(newLVD, rev(newLVD)), border = F, col = "#B1A7926e")
lines(LVDpred$fit ~ newLVD, lwd = 4, col = "white")
text(rec ~ L_VD, labels = Species, data = traits_sp,
     cex = ifelse(traits_sp$Species %in% c("HCR", "TWA"), 1.2, 0.8),
     col = ifelse(traits_sp$Species %in% c("HCR", "TWA"), "black", "grey60"))
mtext(side = 3, adj = 0, line = -1.5, text = " d)", cex = 1.5)
dev.off()


