#########################
#     Trait-Demog       #
# JSDMs with Phylo Dist #
#########################

# run lines 1-110 in Analysis.R first

library(boral)

pd_mat <- cophenetic(tree$scenario.3)

# compare lv correlation structures

null <- boral(y = apply(traits_sp[,25:33],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4, tPC5, tPC6, tPC7, tPC8), 
              family = "normal", calc.ics = T, 
              mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
              lv.control = list(num.lv = 2, type = "independent", distmat = NULL)
)

phylo.exp <- boral(y = apply(traits_sp[,25:33],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4, tPC5, tPC6, tPC7, tPC8),
                   family = "normal", calc.ics = T,
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "exponential", distmat = pd_mat)
)

phylo.sqexp <- boral(y = apply(traits_sp[,25:33],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4, tPC5, tPC6, tPC7, tPC8), 
                     family = "normal", 
                     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                     lv.control = list(num.lv = 2, type = "squared.exponential", distmat = pd_mat)
)

phylo.powexp <- boral(y = apply(traits_sp[,25:33],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4, tPC5, tPC6, tPC7, tPC8), 
                      family = "normal", 
                      mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                      lv.control = list(num.lv = 2, type = "powered.exponential", distmat = pd_mat)
)

phylo.spher <- boral(y = apply(traits_sp[,25:33],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4, tPC5, tPC6, tPC7, tPC8), 
                     family = "normal", 
                     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                     lv.control = list(num.lv = 2, type = "spherical", distmat = pd_mat)
)


offset <- c(-0.3, -0.15, 0, 0.15, 0.3)
library(viridis)
cols <- viridis(5)
colt <- adjustcolor(cols, alpha.f = 0.25)
nt <- 8
nd <- 9
title <- c(" a) Trait PC1",
           " b) Trait PC2",
           " c) Trait PC3",
           " d) Trait PC4",
           " e) Trait PC5",
           " f) Trait PC6",
           " g) Trait PC7",
           " h) Trait PC8")

#pdf("D:/Dropbox/Functional Traits Project/Figures/Boral coef plot Dec21.pdf", height = 12, width = 10)
par(mfrow=c(4,2), mar=c(4,4,2,2), oma=c(1,7,1,1))
for(i in 1:8){
    plot(x = 1:nd, y = 1:nd, type = "n", cex.axis = 1.5,
         xlim = c(-4,4), ylim = c(0.5,nd+0.5), xlab = "", yaxt = "n", ylab = "")
    abline(v = 0, lty = 2)
    for(j in 1:5){
        if(j==1) mod <- null
        if(j==2) mod <- phylo.exp
        if(j==3) mod <- phylo.sqexp
        if(j==4) mod <- phylo.powexp
        if(j==5) mod <- phylo.spher
        colsig <- rep(colt[j], nd)
        sig.index <- which(mod$hpdintervals$X.coefs[,i,"lower"] > 0 | mod$hpdintervals$X.coefs[,i,"upper"] < 0)
        colsig[sig.index] <- cols[j]
        ##colsig <- which(sig==1)
        points(c(1:nd) + offset[j] ~ mod$X.coefs.median[,i], pch = 16, cex = 2, col = colsig)
        
        for(k in 1:nd){
            arrows(mod$hpdintervals$X.coefs[k,i,"lower"], k + offset[j], mod$hpdintervals$X.coefs[k,i,"upper"], k + offset[j], 
                   length = 0, col = colsig[k])
        }
    }
    if(i %in% c(1,3,5,7)) axis(side = 2, at = 1:nd, labels = names(traits_sp)[25:33], las = 1, cex.axis = 1.5)
    mtext(side = 3, adj = 0, cex = 1.15, text = title[i], line = 1.5)
}
legend('topright', bg = "white", pch = 16, col = cols, pt.cex = 2, legend = c(
    "No phylogenetic structure",
    "Exponential",
    "Squared exponential",
    "Powered exponential",
    "Spherical")
)
mtext(side = 1, outer = T, "Effect size", line = -1, cex = 1.5)
dev.off()


####################
# Trait-Swamp JSDM #
####################

# read in data
trees <- read.csv("./raw_data/plot trees (29 spp).csv")
plots <- read.csv("./raw_data/plot.csv")

# match swamp_area data
trees$swamp_area <- plots$swamp_area[match(trees$plot, plots$plot)]/400

# aggregate to tree level to prevent double counting multi-stemmed trees
library(stringr)
trees$UID <- toupper(trees$specimen_no.)
trees$UID <- ifelse(str_ends(trees$UID, "[[A-Z]]"),
                    str_sub(trees$UID, 1, -2),
                    trees$UID)

# only use 2018 data 
trees <- trees[-which(trees$DBH_2018==""|trees$DBH_2018=="cnf"|trees$DBH_2018=="dead"|trees$DBH_2018=="missed"),]
trees$DBH_2018 <- as.numeric(trees$DBH_2018)

trees2018 <- aggregate(DBH_2018 ~ plot + swamp_area + species, 
                       data = trees,
                       FUN = function(x) {sqrt(sum(x^2, na.rm=T))})
trees2018$ba <- pi*(trees2018$DBH_2018/2)^2

library(reshape2)
plottrees <- dcast(trees2018, plot + swamp_area ~ species, fun.aggregate = sum, value.var = "ba", na.rm = T)
Y <- plottrees[, 3:length(plottrees)]
X <- data.frame(swamp_area = plottrees$swamp_area)

# reorder traits so that they match the order of Y's columns
abbrev <- function(x){
    y <- gsub("\\.", " ", x)
    toupper(paste0(
        substr(lapply(strsplit(y, " "), "[[", 1), 0, 1),
        substr(sapply(strsplit(y, " "), "[[", 2), 0, 2)
    ))}
TR <- traits_datonly[match(abbrev(names(Y)), traits_sp$Species),]
TR.PC <- PCA.traits$CA$u[match(abbrev(names(Y)), traits_sp$Species), 1:8]


dim(Y); dim(X); dim(TR)
dim(antTraits$abun); dim(antTraits$env); dim(antTraits$traits)


library(gllvm)
# null model without traits or env
jsdm.null <- gllvm(y = Y, family = "tweedie", num.lv = 2)
plot(jsdm.null)

# model with only env (swamp)
jsdm.swamp <- gllvm(y = Y, X = X, 
                    formula = ~ swamp_area,
                    family = "tweedie", num.lv = 2)
coefplot(jsdm.swamp)

# model with traits. too many traits
jsdm.tr <- gllvm(y = Y, X = X, TR = TR, 
                 family = "tweedie", num.lv = 2)
coefplot(jsdm.tr)

# model with trait PCs
jsdm.trpc <- gllvm(y = Y, X = X, TR = TR.PC, 
                   family = "tweedie", num.lv = 2)
coefplot(jsdm.trpc)

jsdm.trpc.1.1 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,1]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.2 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,2]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.3 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,3]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.4 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,4]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.5 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,5]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.6 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,6]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.7 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,7]), 
                       family = "tweedie", num.lv = 2)
jsdm.trpc.1.8 <- gllvm(y = Y, X = X, TR = data.frame(TR.PC[,8]), 
                       family = "tweedie", num.lv = 2)
which.min(
    AICc.gllvm(jsdm.null,
               jsdm.swamp,
               jsdm.tr,
               jsdm.trpc,
               jsdm.trpc.1.1,
               jsdm.trpc.1.2,
               jsdm.trpc.1.3,
               jsdm.trpc.1.4,
               jsdm.trpc.1.5,
               jsdm.trpc.1.6,
               jsdm.trpc.1.7,
               jsdm.trpc.1.8
    ))
AIC(jsdm.null, jsdm.trpc)



# Clustering traits?
traits.clust <- hclust(dist(traits_datonly))
plot(traits.clust, labels = traits_sp$Species)
