library(cluster)


# start by looking at "high" correlation in corMvn simulation

load("sim50Mvn.sol.highCor.RData", verbose = T)

# look at one individual sim to see how it reacts to diff stats

curSol <- sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]]
curSol <- sim50Mvn.sol.highCor$sol[[43]]$allBcSol[[1]]
length(sim50Mvn.sol.highCor$sol[[42]]$allBcSol[[1]])



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
               list(clusters = list(postProcess(aSol$allBcSol[[1]], ...)))
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

debugonce(postSubSample.tight)
postSubSample.tight(curSol)


debugonce(postSubSample.pca)
postSubSample.pca(curSol)


boxplot(list(orig = computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp),
             new50 = computePerf(cor50Mvn.truth, kmeans.50, amrs.hp),
             new75 = computePerf(cor50Mvn.truth, kmeans.75, amrs.hp),
             mean = computePerf(cor50Mvn.truth, kmeans.mean, amrs.hp),
             tight = computePerf(cor50Mvn.truth,tightMVN, amrs.hp),
             tightMed = computePerf(cor50Mvn.truth, tightMVN.med, amrs.hp),
             tight2 = computePerf(cor50Mvn.truth,tightMVN2, amrs.hp),
             # pca = computePerf(cor50Mvn.truth, pca1, amrs.hp),
             pcaOnly = computePerf(cor50Mvn.truth, pcaOnly, amrs.hp)
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

gasPSS.50.2 <- postSubSample.median.kmeans2(bc30.fdr30p90Z)
meow <- postSubSample.tight(bc30.fdr30p90Z, 0.5)

plotHeatmap(fdr30p90 , gasPSS.50.2, T, T)
plotHeatmap(fdr30p90Z , meow, T, T)
plotHeatmap(fdr30p90 , meow, T)

plotClusterExpression(fdr30p90Z , meow)

fdr30p90[meow$rowIdx, ]


