library(cluster)
pca32 <- postSubSample.pca(condOpt.29.33$sol[[4]])
pca32$dDat$rotation[,"norm1"] <- pca32$dDat$rotation[,"PC1"] * pca32$dDat$sdev[1]
pca32$dDat$rotation[,"norm2"] <- pca32$dDat$rotation[,"PC2"] * pca32$dDat$sdev[2]
ggplot(pca32$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()
library(gridExtra)


# start by looking at "high" correlation in corMvn simulation

load("sim50Mvn.sol.highCor.RData", verbose = T)
load("~/bc.RData")

# look at one individual sim to see how it reacts to diff stats

curSol <- sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]]
mat <- sim50Mvn.sol.highCor$data[[42]]$mat
curSol <- sim50Mvn.sol.highCor$sol[[43]]$allBcSol[[1]]
length(sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]])

tmpSol <- postSubSample.pca(curSol)

simPlot <- corFigure(mat, tmpSol, T)
simPlot <- simPlot + theme(legend.position = c(0.925,0.85), legend.title = element_blank())

simPlot

simPlot + guides(fill = guide_legend(label.position = "left")) + opts(legend.position = c(1,1))

flyPlot <- corFigure(flyWorm, tmpPca, F)
flyPlot <- flyPlot + theme(legend.position = "none")

grid.arrange(simPlot, flyPlot, ncol = 2)
dev.print(pdf, "~/Dropbox/biclustering/ismb/extAbstract/density.pdf")

dev.print(pdf, "~/Dropbox/biclustering/ismb/extAbstract/hist.pdf")


ct <- cor(t(mat[1:300, 1:30]))
plot(density(ct[upper.tri(ct)]))

corFigure <- function(data, sol, truth = FALSE)
{
    corClust <- cor(t(data[sol$rowIdx, sol$colIdx]))
    corAll <- cbind(corClust[upper.tri(corClust)], "Bicluster")

    ranRow <- sample.int(n = nrow(data), size = length(sol$rowIdx))
    ranRowCor <- cor(t(data[ranRow, sol$colIdx]))
    corAll <- rbind(corAll, cbind(ranRowCor[upper.tri(ranRowCor)], "Random rows"))


    ranRow <- sample.int(n = nrow(data), size = length(sol$rowIdx))
    ranCol <- sample.int(n = ncol(data), size = length(sol$colIdx))
    bothRanCor <- cor(t(data[ranRow, ranCol]))
    corAll <- rbind(corAll, cbind(bothRanCor[upper.tri(bothRanCor)], "Both random"))

    allCond <- cor(t(data[sol$rowIdx,]))
    corAll <- rbind(corAll, cbind(allCond[upper.tri(allCond)], "All conditions"))

    if (truth)
    {
        corTruth <- cor(t(data[1:300, 1:30]))
        corAll <- rbind(corAll, cbind(corTruth[upper.tri(corTruth)], "Truth"))
    }

    corAll <- as.data.frame(corAll, stringsAsFactors = F)
    corAll[,1] <- as.numeric(corAll[,1])
    colnames(corAll) <- c("Correlation", "Conditions")
    
    # ggplot(corAll, aes(x = Conditions, y = abs(Correlation), colour = Conditions)) + 
    # ggplot(corAll, aes(x = Conditions, y = abs(Correlation), colour = Conditions)) + 
    ggplot(corAll, aes(x = abs(Correlation), y = ..density.., colour = Conditions)) + 
        # geom_density(aes(fill = Conditions), alpha = 0.5) + 
        geom_histogram(aes(fill = Conditions), position = "dodge") + 
        # geom_boxplot(aes(fill = Conditions), alpha = 0.5) + 
        scale_fill_manual(values=cbbPalette) + 
        scale_colour_manual(values=cbbPalette)
}




sol.50 <- postSubSample.median(curSol)
sol.75 <- postSubSample.median(curSol, 0.75)


plot(sol$)
par(mfrow = c(2,1))
hist(sol.50$ab)
hist(sol.75$ab)

debug(postSubSample.median)

# using gap statistic for clustering

medKmeans.50 <- postSubSample.median.kmeans(curSol)
medKmeans.75 <- postSubSample.median.kmeans(curSol, 0.75, 0.75)

testPostMethod <- function(simData, postProcess, ...)
{
    lapply(simData$sol, function(aSol)
           {
               res <- NA
               tryCatch(
                        res <- postProcess(aSol$allBcSol[[1]], ...), 
                        error = function(e) cat("*****Err: \n\t", 
                                                as.character(e), "\n")
                        )
               # res <- postProcess(aSol$allBcSol[[1]], ...)


               list(clusters = list(res))
           })
}

kmeans.50 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.median.kmeans2)
kmeans.75 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.median.kmeans2, 0.75, 0.75)
kmeans.mean <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.mean.kmeans)
tightMVN <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.tight)
tightMVN2 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.tight)
tightMVN.med <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.tight)

pca1 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca)
pca2 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca)
pcaOnly <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca)

pca90 <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca, 0.90)

pcaCov <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca)

pcaScaled <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.pca)

tPca <- testPostMethod(sim50Mvn.sol.highCor, postSubSample.tightPCA)
tmp <- which(sapply(tPca, function(x) is.na(x$clusters)))
tPca[tmp] <- NULL

tryCatch({
    stop("meow")}, error = function (e) cat(as.character(e), "\nfuuuu")
    )

debugonce(postSubSample.tight)
postSubSample.tight(curSol)
tmp <- postSubSample.tightPCA(curSol)

tmp <- postSubSample.pca(curSol)
tmp <- postSubSample.pca(curSol, abThresh = 0.7)
tmp <- postSubSample.pca(curSol, abThresh = 0.9)


debugonce(postSubSample.pca)
postSubSample.pca(curSol)


boxplot(list(orig = computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp),
             # new50 = computePerf(cor50Mvn.truth, kmeans.50, amrs.hp),
             # new75 = computePerf(cor50Mvn.truth, kmeans.75, amrs.hp),
             mean = computePerf(cor50Mvn.truth, kmeans.mean, amrs.hp),
             tight = computePerf(cor50Mvn.truth,tightMVN, amrs.hp),
             tightMed = computePerf(cor50Mvn.truth, tightMVN.med, amrs.hp),
             pcaOnly = computePerf(cor50Mvn.truth, pcaOnly, amrs.hp),
             pca90 = computePerf(cor50Mvn.truth, pca90, amrs.hp),
             pcaScaled = computePerf(cor50Mvn.truth, pcaScaled, amrs.hp)
             # pcaCov = computePerf(cor50Mvn.truth, pcaCov, amrs.hp)
             # tPca = computePerf(cor50Mvn.truth, tPca, amrs.hp)
             #tight2 = computePerf(cor50Mvn.truth,tightMVN2, amrs.hp),
             # pca = computePerf(cor50Mvn.truth, pca1, amrs.hp),
             ), ylim = c(0,1))

debugonce(computePerf)

computePerf(cor50Mvn.truth, kmeans.50, amrs.hp)

library("fpc")
prediction.strength(as.matrix(sol.50$ab), Gmin = 2, Gmax = 4, clustermethod=kmeansCBI,
                    classification="centroid")

as.matrix(sol.50$ab)

## look at  const mean

cmPSS.50 <- testPostMethod(constMean.sol, postSubSample.median.kmeans)
cmPSS.75 <- testPostMethod(constMean.sol, postSubSample.median.kmeans, 0.75, 0.75)
cmPSS.mean <- testPostMethod(constMean.sol, postSubSample.mean.kmeans)


curSol <- constMean.sol$sol[[1]]$allBcSol[[1]]
postSubSample.median.kmeans(curSol)


boxplot(list(orig = computePerf(constMean.truth, constMean.sol$sol, amrs.hp),
             new50 = computePerf(constMean.truth, cmPSS.50, amrs.hp),
             new75 = computePerf(constMean.truth, cmPSS.75, amrs.hp),
             mean = computePerf(constMean.truth, cmPSS.mean, amrs.hp)), ylim = c(0,1))


# pca method...
debug(postSubSample.pca)
postSubSample.pca(curSol)


################################################################################
# look at fly-worm data
################################################################################

debugonce(postSubSample.median.kmeans)
fwPSS.50 <- postSubSample.median.kmeans(fwBC30)

debugonce(postSubSample.median.kmeans2)
fwPSS.50.2 <- postSubSample.median.kmeans2(fwBC30)

ggPlotExpression(flyWorm[fwPSS.50.2$rowIdx, fwPSS.50.2$colIdx])

heatmap(flyWorm[fwPSS.50.2$rowIdx, fwPSS.50.2$colIdx], Colv = NA, col = redgreen(256))

fwPSS.75 <- postSubSample.median.kmeans(fwBC30, 0.75, 0.75)
fwPSS.mean <- postSubSample.mean.kmeans(fwBC30)

debug(postSubSample.pca)
fwPSS.pca <- postSubSample.pca(fwBC30)

undebug(postSubSample.pca)

# tmp
rot <- postSubSample.pca(fwBC30)
rot$abDat$clust <- "bg"
rot$abDat$clust[rot$cluster$rowIdx] <- "bc"
ggplot(rot$abDat, aes(x = PC1, y = PC2, colour = clust)) + geom_point()

# tight.clust is doing something really stupid here... take a look at plot
debug(postSubSample.tightPCA)

################################################################################
# figures to show Prof. Huang
################################################################################
tmp <- postSubSample.tightPCA(fwBC30)

ggplot(tmp$abDat, aes(x = PC1, y = PC2, colour = cluster)) + geom_point()
ggplot(tmp$abDat, aes(x = PC1, y = zeros, colour = cluster)) + geom_point()
ggplot(tmp$abDat, aes(x = PC1, y = cov, colour = cluster)) + geom_point()
ggplot(tmp$abDat, aes(x = median, y = cov, colour = cluster)) + geom_point()
ggplot(tmp$abDat, aes(x = PC1, colour = cluster)) + geom_density()
ggplot(tmp$abDat, aes(x = PC1, colour = cluster)) + geom_density(aes(fill = cluster), alpha = 0.5)
ggplot(tmp$abDat, aes(x = median, colour = cluster)) + geom_density(aes(fill = cluster), alpha = 0.5)

tmpPca <- postSubSample.pca(fwBC30)

# figure for abstract
corFW <- cor(t(flyWorm[tmpPca$rowIdx, tmpPca$colIdx]))
corGG <- cbind(corFW[upper.tri(corFW)], "bc")

set.seed(42)

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette 
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


# To use for fills, add
scale_fill_manual(values=cbbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbbPalette)

# get random genes
ranRow <- sample.int(n = nrow(flyWorm), size = length(tmpPca$rowIdx))
cor2 <- cor(t(flyWorm[ranRow, tmpPca$colIdx]))
corGG <- rbind(corGG, cbind(cor2[upper.tri(cor2)], "ranRow"))

# # get random conditions
# ranCol <- sample.int(n = ncol(flyWorm), size = length(tmpPca$colIdx))
# cor3 <- cor(t(flyWorm[tmpPca$rowIdx, ranCol]))
# corGG <- rbind(corGG, cbind(cor3[upper.tri(cor3)], "ranCol"))

# get random genes and conditions
ranRow <- sample.int(n = nrow(flyWorm), size = length(tmpPca$rowIdx))
ranCol <- sample.int(n = ncol(flyWorm), size = length(tmpPca$colIdx))
cor4 <- cor(t(flyWorm[ranRow, ranCol]))
corGG <- rbind(corGG, cbind(cor4[upper.tri(cor4)], "bothRan"))

cor5 <- cor(t(flyWorm[tmpPca$rowIdx,]))
corGG <- rbind(corGG, cbind(cor5[upper.tri(cor5)], "allCond"))

corGG <- as.data.frame(corGG, stringsAsFactors = FALSE)
corGG[,1] <- as.numeric(corGG[,1])
colnames(corGG) <- c("Correlation", "Type")

ggplot(corGG, aes(x = abs(Correlation), colour = Type)) + geom_density(aes(fill = Type), alpha = 0.5) + 
    scale_fill_manual(values=cbbPalette) + 
    scale_colour_manual(values=cbbPalette)

ggsave("../ismb/extAbstract/cor.pdf")
# Saving 8.02 x 3.44 in image

ggplot(corGG, aes(x = abs(Correlation), colour = Type)) + geom_freqpoly(aes(fill = Type), alpha = 0.5)

ggplot(tmpPca$abDat, aes(x = PC1, y = PC2, colour = cluster)) + geom_point()
ggplot(tmpPca$abDat, aes(x = PC1, y = zeros, colour = cluster)) + geom_point()
ggplot(tmpPca$abDat, aes(x = PC1, y = cov, colour = cluster)) + geom_point()
ggplot(tmpPca$abDat, aes(x = median, y = cov, colour = cluster)) + geom_point()
ggplot(tmpPca$abDat, aes(x = PC1, colour = cluster)) + geom_density()
ggplot(tmpPca$abDat, aes(x = median, colour = cluster)) + geom_density()


elust1 <- flybaseNames[which(tmp$abDat$cluster == 1)]
clust2 <- flybaseNames[which(tmp$abDat$cluster == 2)]
clust3 <- flybaseNames[which(tmp$abDat$cluster == 3)]
clust4 <- flybaseNames[which(tmp$abDat$cluster == 4)]

clustUnion <- with(tmp$abDat, which(cluster == 2 | cluster == 3 | cluster == 4))

clust1.worm <- wormNames[which(tmp$abDat$cluster == 1)]
clust2.worm <- wormNames[which(tmp$abDat$cluster == 2)]
clust3.worm <- wormNames[which(tmp$abDat$cluster == 3)]
clust4.worm <- wormNames[which(tmp$abDat$cluster == 4)]

ggPlotExpression(flyWorm[which(tmp$abDat$cluster == 1), tmp$cluster$colIdx])
ggPlotExpression(flyWorm[which(tmp$abDat$cluster == 2), tmp$cluster$colIdx])
ggPlotExpression(flyWorm[which(tmp$abDat$cluster == 3), tmp$cluster$colIdx])
ggPlotExpression(flyWorm[which(tmp$abDat$cluster == 4), tmp$cluster$colIdx])

plotHeatmap2(flyWorm[which(tmp$abDat$cluster == 1), tmp$cluster$colIdx])
plotHeatmap2(flyWorm[which(tmp$abDat$cluster == 2), tmp$cluster$colIdx])
plotHeatmap2(flyWorm[which(tmp$abDat$cluster == 3), tmp$cluster$colIdx])
plotHeatmap2(flyWorm[which(tmp$abDat$cluster == 4), tmp$cluster$colIdx])

# for pcaOnly results
plotHeatmap2(flyWorm[tmpPca$cluster$rowIdx, tmpPca$cluster$colIdx], cols = T)

write(clust1, "../tmp/clust1.txt")
write(clust2, "../tmp/clust2.txt")
write(clust3, "../tmp/clust3.txt")
write(clust4, "../tmp/clust4.txt")
write(flybaseNames, "../tmp/bg.txt")


write(wormNames[tmpPca$cluster$rowIdx], "../tmp/pcaOnly.worm.txt")

write(clust1.worm, "../tmp/clust1.worm.txt")
write(clust2.worm, "../tmp/clust2.worm.txt")
write(clust3.worm, "../tmp/clust3.worm.txt")
write(clust4.worm, "../tmp/clust4.worm.txt")
write(wormNames, "../tmp/bgWorm.txt")


# this clusters using only the first PC
tmp <- postSubSample.pca(fwBC30, abThresh = 0.89)

tmp <- postSubSample.pca(fwBC30)

tmp$abDat$clust <- "bg"
tmp$abDat$clust[tmp$cluster$rowIdx] <- "bc"
ggplot(tmp$abDat, aes(x = PC1, y = PC2, colour = clust)) + geom_point()
ggplot(tmp$abDat, aes(x = PC1, y = cov, colour = clust)) + geom_point()
ggplot(tmp$abDat, aes(x = PC1)) + geom_density(aes(fill = clust), alpha = 0.5)
# median and PC1 correlate well
ggplot(tmp$abDat, aes(x = PC1, y = median, colour = as.factor(clust))) + geom_point()
ggplot(tmp$abDat, aes(x = PC1, y = var, colour = as.factor(clust))) + geom_point()
ggplot(tmp$abDat, aes(x = median, y = var, colour = as.factor(clust))) + geom_point()

# regularizing more doesn't seem to change the picture much
tmp2 <- postSubSample.pca(fw30.reg3)
tmp2$abDat$clust <- "bg"
tmp2$abDat$clust[tmp2$cluster$rowIdx] <- "bc"
ggplot(tmp2$abDat, aes(x = PC1, y = PC2, colour = clust)) + geom_point()
# median and PC1 correlate well
ggplot(tmp2$abDat, aes(x = PC1, y = median, colour = as.factor(clust))) + geom_point()
ggplot(tmp2$abDat, aes(x = PC1, y = var, colour = as.factor(clust))) + geom_point()
ggplot(tmp2$abDat, aes(x = median, y = var, colour = as.factor(clust))) + geom_point()

tReg <- postSubSample.tight(fwBC30)
tReg$abDat$clust <- "bg"
tReg$abDat$clust[tReg$cluster$rowIdx] <- "bc"
ggplot(tReg$abDat, aes(x = med, y = mean, colour = clust)) + geom_point()


################################################################################
## end figures
################################################################################

ggPlotParSolution(fwBC30)
ggPlotParSolution(sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]])

# big problem -- not really separable in 1D, try different regularization

fw30.HL <- bcSubSamplePar(flyWorm, 100, 30, 0.6, 
                          clustOptions = list(penalty = "HL"))

fw30.SCAD <- bcSubSamplePar(flyWorm, 100, 30, 0.6, 
                          clustOptions = list(penalty = "SCAD"))

save(fw30.HL, fw30.SCAD, fwBC30, file = "fwPenalty.RData")


fw30.SOFT <- bcSubSamplePar(flyWorm, 100, 30, 0.6, 
                            clustOptions = list(penalty = "SOFT"))

load("../bcSol/fwPenalty.RData", verbose = T)

debugonce(postSubSample.median.kmeans)
postSubSample.median.kmeans(fw30.HL)

debugonce(postSubSample.mean.kmeans)
postSubSample.mean.kmeans(fw30.HL)

postSubSample.median.kmeans(fw30.SCAD)

# this still didn't work out too well... try more regularization

# getting "optimal" regularization w/ default params
fw30 <- bcSubSamplePar(flyWorm, 10, 30, 0.6)
fw30.lam <- t(sapply(fw30, function(x) x$sccaLam))

# getting "optimal" regularization w/ default params
fw30.reg1 <- bcSubSamplePar(flyWorm, 100, 30, 0.6,
                            clustOptions = list(lamx = seq(0.3, 1, length.out = 4)))
t(sapply(fw30.reg1, function(x) x$sccaLam))

fw30.reg2 <- bcSubSamplePar(flyWorm, 100, 30, 0.6,
                            clustOptions = list(lamx = seq(0.1, 0.8, length.out = 4)))
t(sapply(fw30.reg2, function(x) x$sccaLam))

fw30.reg3 <- bcSubSamplePar(flyWorm, 100, 30, 0.6,
                            clustOptions = list(lamx = seq(0.03, 0.25, length.out = 4)))
t(sapply(fw30.reg3, function(x) x$sccaLam))

fw30.reg4 <- bcSubSamplePar(flyWorm, 100, 30, 0.6,
                            clustOptions = list(lamx = seq(0.01, 0.1, length.out = 4)))
t(sapply(fw30.reg4, function(x) x$sccaLam))

debug(postSubSample.median.kmeans2)

# Doesn't work very well
# reg1 <- postSubSample.median.kmeans2(fw30.reg1)
# reg2 <- postSubSample.median.kmeans2(fw30.reg2)
# reg3 <- postSubSample.median.kmeans2(fw30.reg3)


apply
ggPlotParSolution(fwBC30)

noReg <- postSubSample.tight(fwBC30, 0.6)

debug(postSubSample.tight)
reg1 <- postSubSample.tight(fw30.reg1, 0.6)
reg2 <- postSubSample.tight(fw30.reg2, 0.6)
reg3 <- postSubSample.tight(fw30.reg3, 0.6)
reg4 <- postSubSample.tight(fw30.reg4, 0.6)

reg1.1 <- postSubSample.pca(fw30.reg1)


plotHeatmap(flyWorm, reg1)
plotHeatmap(flyWorm, reg2)
plotHeatmap(flyWorm, reg3)
plotHeatmap(flyWorm, noReg)

vScore(noReg$rowIdx, reg1$rowIdx)
vScore(reg1$rowIdx, reg2$rowIdx)
vScore(reg1$rowIdx, reg3$rowIdx)
vScore(reg2$rowIdx, reg3$rowIdx)
vScore(reg3$rowIdx, reg4$rowIdx)

save(fw30.reg1, fw30.reg2, fw30.reg3, fw30.reg4, file = "fwReg.RData")

load("../bcSol/fwReg.RData", verbose = T)

################################################################################
# look at gene data
################################################################################
load("../bcSol/bc30.fdr30p90Z.RData", verbose = T)

plotHeatmap <- function(data, clust, rows = T, cols = F)
{
    if (rows)
        rows <- NULL
    else
        rows <- NA
    if (cols)
        cols <- NULL
    else
        cols <- NA

    heatmap(data[clust$rowIdx, clust$colIdx], Rowv = rows, Colv = cols, col = redgreen(256))
}

plotHeatmap2 <- function(data, rows = T, cols = F)
{
    if (rows)
        rows <- NULL
    else
        rows <- NA
    if (cols)
        cols <- NULL
    else
        cols <- NA

    heatmap(data, Rowv = rows, Colv = cols, col = redgreen(256))
}

gasPSS.50.2 <- postSubSample.median.kmeans2(bc30.fdr30p90Z)
meow <- postSubSample.tight(bc30.fdr30p90Z, 0.5)

plotHeatmap(fdr30p90 , gasPSS.50.2, T, T)
plotHeatmap(fdr30p90Z , meow, T, T)
plotHeatmap(fdr30p90 , meow, T)

plotClusterExpression(fdr30p90Z , meow)

fdr30p90[meow$rowIdx, ]


# driver for condition seleciton
curMat <- sim50Mvn.sol.highCor$data[[40]]$mat
curMatZ <- t(scale(t(curMat)))

debug(optConditionSize)
condOpt.29.33 <- optConditionSize(curMat, 29, 33)
condOpt.29.33.Z <- optConditionSize(curMatZ, 29, 33)
save(curMat, curMatZ, condOpt.29.33, condOpt.29.33.Z, 
     file = "condOpt.RData")

### Look at all diff solutions & investigate why cond vector is small
load("~/Dropbox/biclustering/bcSol/condOpt.RData", verbose = T)
ggPlotParSolution(condOpt.29.33$sol[[1]])
ggPlotParSolution(condOpt.29.33$sol[[2]])
ggPlotParSolution(condOpt.29.33$sol[[3]])
ggPlotParSolution(condOpt.29.33$sol[[4]])

apply(getD(condOpt.29.33$sol[[1]]), 1, mean)

pca30 <- postSubSample.pca(condOpt.29.33$sol[[2]])
pca30$dDat$rotation[,"norm1"] <- pca30$dDat$rotation[,"PC1"] * pca30$dDat$sdev[1]
pca30$dDat$rotation[,"norm2"] <- pca30$dDat$rotation[,"PC2"] * pca30$dDat$sdev[2]
ggplot(pca30$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()

pca31 <- postSubSample.pca(condOpt.29.33$sol[[3]])
pca31$dDat$rotation[,"norm1"] <- pca31$dDat$rotation[,"PC1"] * pca31$dDat$sdev[1]
pca31$dDat$rotation[,"norm2"] <- pca31$dDat$rotation[,"PC2"] * pca31$dDat$sdev[2]
ggplot(pca31$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()

pca32 <- postSubSample.pca(condOpt.29.33$sol[[4]])
pca32$dDat$rotation[,"norm1"] <- pca32$dDat$rotation[,"PC1"] * pca32$dDat$sdev[1]
pca32$dDat$rotation[,"norm2"] <- pca32$dDat$rotation[,"PC2"] * pca32$dDat$sdev[2]
ggplot(pca32$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()


### Look at Z transformed
ggPlotParSolution(condOpt.29.33.Z$sol[[1]])
ggPlotParSolution(condOpt.29.33.Z$sol[[2]])
ggPlotParSolution(condOpt.29.33.Z$sol[[3]])
ggPlotParSolution(condOpt.29.33.Z$sol[[4]])

pca31 <- postSubSample.pca(condOpt.29.33.Z$sol[[2]])
pca31$dDat$rotation[,"norm1"] <- pca31$dDat$rotation[,"PC1"] * pca31$dDat$sdev[1]
pca31$dDat$rotation[,"norm2"] <- pca31$dDat$rotation[,"PC2"] * pca31$dDat$sdev[2]
ggplot(pca31$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()

pca33 <- postSubSample.pca(condOpt.29.33$sol[[4]])
pca33$dDat$rotation[,"norm1"] <- pca33$dDat$rotation[,"PC1"] * pca33$dDat$sdev[1]
pca33$dDat$rotation[,"norm2"] <- pca33$dDat$rotation[,"PC2"] * pca33$dDat$sdev[2]
ggplot(pca33$dDat$rotation, aes(x = norm1, y = norm2, colour = cluster)) + geom_point()
ggPlotParSolution(condOpt.29.33$sol[[4]])
###

ggPlotExpression(curMatZ[1:600, ])

fullSol <- biclusteringPar(curMat, lam.lwr = 20, lam = 30)
psFullSol <- postSubSample.pca(fullSol)

ssSol <- bcSubSamplePar(curMat, lam.lwr = 20, lam = 30)
save(ssSol, curMat, psSol, fullSol, file = "problemSim.RData")

load("~/Dropbox/biclustering/bcSol/problemSim.RData", verbose= T)

# These are the sims w/ problems
ggPlotParSolution(fullSol)
debugonce(ggPlotParSolution)

ggPlotParSolution(sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]])

psSol <- postSubSample.pca(ssSol)

ssSolZ <- bcSubSamplePar(curMatZ, lam.lwr = 20, lam = 30)

