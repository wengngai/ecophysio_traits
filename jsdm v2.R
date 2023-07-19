#########################
#     Trait-Demog       #
# JSDMs with Phylo Dist #
#########################

# perform JSDM on demographic parameters (as multivariate data that are internally correlated)
# run lines 1-110 in Analysis.R first
sourcePartial <- function(fn, skip=0, n=-1) {
    lines <- scan(fn, what=character(), sep="\n", skip=skip, n=n, quiet=TRUE)
    tc <- textConnection(lines)
    source(tc)
    close(tc)
}
sourcePartial("analysis.R", 0, 110)

library(boral)

pd_mat <- cophenetic.phylo(tree$scenario.3)

# compare lv (phylogenetic) correlation structures

null <- boral(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
              family = "normal", calc.ics = T, save.model = T,
              mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
              lv.control = list(num.lv = 2, type = "independent", distmat = NULL)
)

phylo.exp <- boral(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4),
                   family = "normal", calc.ics = T, save.model = T,
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "exponential", distmat = pd_mat)
)

phylo.sqexp <- boral(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                     family = "normal", save.model = T, 
                     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                     lv.control = list(num.lv = 2, type = "squared.exponential", distmat = pd_mat)
)

phylo.powexp <- boral(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                      family = "normal", save.model = T, 
                      mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                      lv.control = list(num.lv = 2, type = "powered.exponential", distmat = pd_mat)
)

phylo.spher <- boral(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                     family = "normal", save.model = T, 
                     mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                     lv.control = list(num.lv = 2, type = "spherical", distmat = pd_mat)
)

# create a list to store 90% HPD intervals
hpd90 <- list()
hpd90[[1]] <- get.hpdintervals(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                              fit.mcmc = get.mcmcsamples(null), 
                              lv.control = list(num.lv = 2, type = "independent", distmat = NULL), prob = 0.9)
hpd90[[2]] <- get.hpdintervals(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                    fit.mcmc = get.mcmcsamples(phylo.exp), 
                    lv.control = list(num.lv = 2, type = "exponential", distmat = pd_mat), prob = 0.9)
hpd90[[3]] <- get.hpdintervals(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                               fit.mcmc = get.mcmcsamples(phylo.sqexp), 
                               lv.control = list(num.lv = 2, type = "squared.exponential", distmat = pd_mat), prob = 0.9)
hpd90[[4]] <- get.hpdintervals(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                               fit.mcmc = get.mcmcsamples(phylo.powexp), 
                               lv.control = list(num.lv = 2, type = "powered.exponential", distmat = pd_mat), prob = 0.9)
hpd90[[5]] <- get.hpdintervals(y = apply(traits_sp[,26:32],2,scale), X = cbind(tPC1, tPC2, tPC3, tPC4), 
                               fit.mcmc = get.mcmcsamples(phylo.spher), 
                               lv.control = list(num.lv = 2, type = "spherical", distmat = pd_mat), prob = 0.9)

offset <- c(-0.34, -0.17, 0, 0.17, 0.34)
cols <- c("#FC7753", "#66D7D1", "#403D58", "#8E8963", "#DBD56E")
colt <- adjustcolor(cols, alpha.f = 0.25)
nt <- 4
nd <- 7
title <- c(" a) Trait PC1",
           " b) Trait PC2",
           " c) Trait PC3",
           " d) Trait PC4")
names(traits_sp[,26:32])
axislabels <- c("Growth: a",
                "Growth: b",
                "Growth: c",
                "Survival: K",
                "Survival: r",
                "Survival: p",
                "Sapling recruitment: BNR")


#pdf("D:/Dropbox/Functional Traits Project/Figures/Boral coef plot Dec21 traits 1-4 v2.pdf", height = 14, width = 12)
par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(1,15,1,1))
for(i in 1:4){
    plot(x = 1:nd, y = 1:nd, type = "n", cex.axis = 1.5,
         xlim = c(-5,5), ylim = c(0.5,nd+0.5), xlab = "", yaxt = "n", ylab = "")
    abline(v = 0, lty = 2)
    for(j in 1:5){
        if(j==1) mod <- null
        if(j==2) mod <- phylo.exp
        if(j==3) mod <- phylo.sqexp
        if(j==4) mod <- phylo.powexp
        if(j==5) mod <- phylo.spher
        colsig <- rep(colt[j], nd)
        sig.index <- which(hpd90[[j]]$X.coefs[,i,"lower"] > 0 | hpd90[[j]]$X.coefs[,i,"upper"] < 0)
        colsig[sig.index] <- cols[j]
        points(c(1:nd) + offset[j] ~ mod$X.coefs.median[,i], pch = 16, cex = 2, col = colsig)
        
        for(k in 1:nd){
            arrows(mod$hpdintervals$X.coefs[k,i,"lower"], k + offset[j], mod$hpdintervals$X.coefs[k,i,"upper"], k + offset[j], 
                   length = 0, col = colsig[k])
            arrows(hpd90[[j]]$X.coefs[k,i,"lower"], k + offset[j], hpd90[[j]]$X.coefs[k,i,"upper"], k + offset[j], 
                   lwd = 3, length = 0, col = colsig[k])
            
        }
    }
    if(i %in% c(1,3,5,7)) axis(side = 2, at = 1:nd, labels = axislabels, las = 1, cex.axis = 1.5)
    mtext(side = 3, adj = 0, cex = 1.15, text = title[i], line = 0.5)
    if(i == 4){
        legend('bottomright', bg = "white", pch = 16, col = cols, pt.cex = 2, legend = c(
            "No phylogenetic structure",
            "Exponential",
            "Squared exponential",
            "Powered exponential",
            "Spherical")
        )
    }
}
mtext(side = 1, outer = T, "Effect size", line = -1, cex = 1.5)
dev.off()


# Note to self: species locations look strange. Try a null boral model with NO X variables (TRAIT PCs) then plot lvsplot
# in this plot, growth-survival tradeoff is evident. maybe without traits that will disappear

phylo.exp.null <- boral(y = apply(traits_sp[,26:32],2,scale),
                   family = "normal", calc.ics = T, save.model = T,
                   mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 50),
                   lv.control = list(num.lv = 2, type = "exponential", distmat = pd_mat)
)


lvvals <- lvsplot(phylo.exp, return.vals = T)
lvvals.null <- lvsplot(phylo.exp.null, return.vals = T)

vec <- lvvals$scaled.lv.coefs
vec.null <- lvvals.null$scaled.lv.coefs
vec.null[,2] <- -1 * vec.null[,2]
rownames(phylo.exp$lv.coefs.median)[8] <- "BSR"

par(mfrow = c(1,2))
plot(1.2*vec, type = "n", )
text(1.4 * lvvals$scaled.lvs, labels = traits_sp$Species, col = "grey")
for(i in 1:nrow(vec)){
    arrows(0, 0, vec[i,1], vec[i,2], col = "red")
    text(vec[i,1] * 1.1, vec[i,2] * 1.1, labels = rownames(phylo.exp$lv.coefs.median)[i], col = "red")
}
plot(1.4*vec.null, type = "n", )
text(1.4 * lvvals.null$scaled.lvs, labels = traits_sp$Species, col = "grey")
for(i in 1:nrow(vec.null)){
    arrows(0, 0, vec.null[i,1], vec.null[i,2], col = "red")
    text(vec.null[i,1] * 1.1, vec.null[i,2] * 1.1, labels = rownames(phylo.exp$lv.coefs.median)[i], col = "red")
}


expcors <- get.residual.cor(phylo.exp, est = "median", prob = 0.9)
library(corrplot)
corrplot(expcors$sig.cor, diag = F, type = "lower")


