set.seed(42)
figBlock <- generateRandomBlock.progRows(50, 20)
# figBlock <- figBlock$sample
corFigBlock <- cor(t(figBlock))

bpColors <- colors(F)[c(35, 50, 128, 655)]
bpColors2 <- colors(F)[c(35, 50, 128, 111, 655)]

hist(corFigBlock[upper.tri(corFigBlock)], 
     prob = TRUE,
     xlab = "Correlation", 
     # main = "Correlation between features in progressive bicluster",
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corProgClust.pdf")

figBlock.iid <- genBlock.iid.N(50, 20, 0, 1)
corFigBlock.iid <- cor(t(figBlock.iid))

hist(corFigBlock.iid[upper.tri(corFigBlock.iid)],
     prob = TRUE,
     xlab = "Correlation", 
     # main = "Correlation between features in indepdent samples",
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corIid.pdf")

set.seed(42)

mvnCovMatrix <- genRanPosDefMat(300, 0.8, 0.95)
plot(mvnCovMatrix[upper.tri(mvnCovMatrix)])

mvnCovMatrix <- genRanPosDefMat(300)
hist(mvnCovMatrix[upper.tri(mvnCovMatrix)],
     prob = TRUE,
     xlab = "Correlation",
     # main = "Sample covariance/correlation in MVN simulation",
     main = "",
     xlim = c(0, 1), col = "lightgrey")
dev.print(pdf, "../img/corMvn5to8.pdf")

mvnCovMatrix <- genRanPosDefMat(300, 0.7, 0.9)
lwr <- seq(0, 1, by = 0.01)

unifVar <- expand.grid(lwr = seq(0, 1, by = 0.01), upr = seq(lwr, 1, by = 0.01))
gridCov <- apply(unifVar, 1, function(row)
                 {
                     if (row[1] >= row[2])
                         return(c(NA, NA))
                     mat <- genRanPosDefMat(300, row[1], row[2])
                     range(mat[upper.tri(mat)])
                 })
gridCov <- t(gridCov)

plot(gridCov[9400:9600,1], col = "blue", pch = 20, lwd = 0.25)
points(gridCov[9400:9600,2], col = "red", pch = 20, lwd = 0.25)

mvnCovMatrix <- genRanPosDefMat(300, 0.8, 0.95)

par(mfrow = c(3,1))

mvnCovMatrix <- genRanPosDefMat(300, .72, 0.94)
range(mvnCovMatrix[upper.tri(mvnCovMatrix)])

hist(mvnCovMatrix[upper.tri(mvnCovMatrix)],
     prob = TRUE,
     xlab = "Correlation",
     # main = "Sample covariance/correlation in MVN simulation U[0.72, 0.94]",
     main = "",
     xlim = c(0, 1), col = "lightgrey")
dev.print(pdf, "../img/corMvn72to94.pdf")

mvnCovMatrix <- genRanPosDefMat(300, .8, 0.95)
range(mvnCovMatrix[upper.tri(mvnCovMatrix)])
hist(mvnCovMatrix[upper.tri(mvnCovMatrix)],
     prob = TRUE,
     xlab = "Correlation",
     main = "Sample covariance/correlation in MVN simulation",
     xlim = c(0, 1))

mvnCovMatrix <- genRanPosDefMat(300)

hist(mvnCovMatrix[upper.tri(mvnCovMatrix)],
     prob = TRUE,
     xlab = "Correlation",
     main = "Sample covariance/correlation in MVN simulation U[0.5, 0.8]",
     xlim = c(0, 1))
dev.print(pdf, "../img/corMvn5to8.pdf")

set.seed(42)
mvnBicluster <- genRandomBlock.mvn(300, 30)
mvnSimDat <- generateNormal(nrow = 1500, ncol = 100, clusterOptions = 
                            list(list(x.start = 1, x.end = 300, y.start = 1, y.end = 30,
                                      data = mvnBicluster)))

ggPlotExpression((mvnSimDat[1:500, 1:100]))
ggsave("../img/corMvnExpression5to8.pdf", width = 21.6, height = 9.91)
levelplot(mvnSimDat[1:500, 1:100], xlab = "Feature" , ylab = "Condition")
dev.print(pdf, "../img/corMvnExpression5to8.pdf")

# boxplots for const mean
load("constMeanComp.RData")
load("constMean.sol.RData")
load("sparseBC.sim1.matrix.sol.RData")

boxplot(
        computePerf(constMean.truth, constMean.sol$sol, amrs.hp),
        computePerf(constMean.truth, plaid.constMean.sol, amrs.hp),
        computePerf(constMean.truth, sparseBC.constMean.sol, amrs.hp),
        computePerf(constMean.truth, sparseBC.sim1.matrix.sol, amrs.hp),
        computePerf(constMean.truth, ssvd.constMean.sol, amrs.hp),
        ylim = c(0, 1), xlab = "Method", 
        ylab = "eAMRS", col = bpColors2)
axis(1, at = 1:5, labels = c("SCCA-BC", "IP", "sparseBC", "matrixBC", "SSVD"))

dev.print(pdf, "../img/constMeanBoxplot.pdf")


mean(sapply(constMean.sol$sol, function(x)
       {
           all.equal(1:15, sort(x$clusters[[1]]$colIdx))
       }))

sort([[3]]$clusters[[1]]$colIdx)

# Boxplots for mvn sim
load("mvn50Comp.RData")
ssvd.sim50mvn.perf <- computePerf(cor50Mvn.truth, ssvd.sim50Mvn, amrs.hp)
plaid.sim50mvn.perf <- computePerf(cor50Mvn.truth, plaid.sim50Mvn, amrs.hp)
simMvn50.perf

boxplot(simMvn50.perf, 
        plaid.sim50mvn.perf, 
        computePerf(cor50Mvn.truth, sparseBC.sim3.normal.sol, amrs.hp),
        ssvd.sim50mvn.perf, 
        ylim = c(0, 1), xlab = "Method",
        ylab = "eAMRS",
        col = bpColors)
axis(1, at = 1:4, labels = c("SCCA-BC", "IP", "sparseBC", "SSVD"))
dev.print(pdf, "../img/mvn5to8Boxplot.pdf")

sBC.sim3a.mat <- computePerf(cor50Mvn.truth, sparseBC.sim3.matrix, amrs.hp)
save(sBC.sim3a.mat, file = "sBC.sim3a.mat.RData")


sBC.sim3b.mat <- computePerf(cor50Mvn.truth, sparseBC.sim3.matrix.highCor.sol, amrs.hp)
save(sBC.sim3b.mat, file = "sBC.sim3b.mat.RData")

load("sBC.sim3a.mat.RData")

boxplot(simMvn50.perf, 
        plaid.sim50mvn.perf, 
        computePerf(cor50Mvn.truth, sparseBC.sim3.normal.sol, amrs.hp),
        sBC.sim3a.mat,
        ssvd.sim50mvn.perf, 
        ylim = c(0, 1), xlab = "Method",
        ylab = "eAMRS",
        col = bpColors2)
axis(1, at = 1:5, labels = c("SCCA-BC", "IP", "sparseBC", "matrixBC", "SSVD"))
dev.print(pdf, "../img/mvn5to8Boxplot.pdf")

load("sim50Mvn.sol.highCor.RData")
load("mvn50comp.highCor.RData")
load("sparseBC.sim3.normal.RData")
load("sparseBC.sim3b.normal.RData")
load("sparseBC.sim3.matrix.RData")
load("sparseBC.sim3.matrix.highCor.RData")
load("sBC.sim3b.mat.RData")

boxplot(
        computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp),
        computePerf(cor50Mvn.truth, plaid.sim50Mvn.highCor, amrs.hp),
        computePerf(cor50Mvn.truth, sparseBC.sim3b.normal.sol, amrs.hp),
        sBC.sim3b.mat,
        computePerf(cor50Mvn.truth, ssvd.sim50Mvn.highCor, amrs.hp),
        ylim = c(0, 1), xlab = "Method", 
        ylab = "eAMRS",
        col = bpColors2)
axis(1, at = 1:5, labels = c("SCCA-BC", "IP", "sparseBC", "matrixBC", "SSVD"))
dev.print(pdf, "../img/mvn72to94Boxplot.pdf")

order(computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp))[25] # idx near the median is 9

nearMedMvn <- sim50Mvn.sol.highCor$sol[[9]]$allBcSol[[1]]
nearMedMvn <- lapply(nearMedMvn, function(x){
                     x$ab <- x$ab[c(1:300, sample(300:1500, 200))]
                     x
                           })
plotParSolution(nearMedMvn, "mvn72to94")

ggPlotParSolution(nearMedMvn, "mvn72to94")

plotParSolution(nearMedMvn, "mvn72to94")

# boxplots for two BC
load("twoBC.sol.RData")
load("twoBC.sol.comp.RData")
load("sparseBC.sim2.RData")
load("sparseBC.sim2.matrix.sol.RData")
truth.twoBC <- list(list(rowIdx = 1:50, colIdx = 1:15), 
                    list(rowIdx = 200:239, colIdx = 26:40))

# this is the "old" version
boxplot(
        computePerf(truth.twoBC, twoBC.sol$sol, amrs.eachClust)[1,],
        computePerf(truth.twoBC, plaid.twoBC.sol, amrs.eachClust)[1,],
        computePerf(truth.twoBC, ssvd.twoBC.sol, amrs.eachClust)[1,],
        ylim = c(0, 1), xlab = "Method", 
        ylab = "eAMRS" )
axis(1, at = 1:3, labels = c("SCCA-BC", "IP", "SSVD"))

# this is the "old" version
boxplot(
        computePerf(truth.twoBC, twoBC.sol$sol, amrs.eachClust)[2,],
        computePerf(truth.twoBC, plaid.twoBC.sol, amrs.eachClust)[2,],
        computePerf(truth.twoBC, ssvd.twoBC.sol, amrs.eachClust)[2,],
        ylim = c(0, 1), xlab = "Method", 
        ylab = "eAMRS" )
axis(1, at = 1:3, labels = c("SCCA-BC", "IP", "SSVD"))
dev.print(pdf, "../img/twoBC2Boxplot.pdf")


par(mfrow = c(1,2))

boxplot(
        computePerf(list(truth.twoBC[[1]]), getClusterX(twoBC.sol$sol, 1), amrs.hp),
        computePerf(list(truth.twoBC[[1]]), getClusterX(plaid.twoBC.sol, 1), amrs.hp),
        computePerf(list(truth.twoBC[[1]]), getClusterX(sparseBC.sim2, 1), amrs.hp),
        computePerf(list(truth.twoBC[[1]]), getClusterX(sparseBC.sim2.matrix.sol, 1), amrs.hp),
        computePerf(list(truth.twoBC[[1]]), getClusterX(ssvd.twoBC.sol, 1), amrs.hp),
        ylim = c(0, 1), xlab = "Method", ylab = "eAMRS",
        col = bpColors2)
axis(1, at = 1:5, labels = c("SCCA-BC", "IP", "sparseBC", "matrixBC", "SSVD"))
dev.print(pdf, "../img/twoBC1Boxplot.pdf")

boxplot(
        computePerf(list(truth.twoBC[[2]]), getClusterX(twoBC.sol$sol, 2), amrs.hp),
        computePerf(list(truth.twoBC[[2]]), getClusterX(plaid.twoBC.sol, 2), amrs.hp),
        computePerf(list(truth.twoBC[[2]]), getClusterX(sparseBC.sim2, 2), amrs.hp),
        computePerf(list(truth.twoBC[[2]]), getClusterX(sparseBC.sim2.matrix.sol, 2), amrs.hp),
        computePerf(list(truth.twoBC[[2]]), getClusterX(ssvd.twoBC.sol, 2), amrs.hp),
        ylim = c(0, 1), xlab = "Method", 
        ylab = "eAMRS", col = bpColors2)
axis(1, at = 1:5, labels = c("SCCA-BC", "IP", "sparseBC", "matrixBC", "SSVD"))
dev.print(pdf, "../img/twoBC2Boxplot.pdf")


order(computePerf(truth.twoBC, twoBC.sol$sol, amrs.eachClust)[2,])[25] # index 34 is the closest to median
plotParSolution(twoBC.sol$sol[[34]]$allBcSol[[2]], "twoBC2")

ggPlotParSolution(twoBC.sol$sol[[34]]$allBcSol[[2]], "twoBC2")

plotParSolution(constMean.sol$sol[[34]]$allBcSol[[1]])

ggPlotParSolution(constMean.sol$sol[[34]]$allBcSol[[1]])


A <- abs(sapply(constMean.sol$sol[[34]]$allBcSol[[1]], function (x) x$ab))
A <- abs(sapply(nearMedMvn, function (x) x$ab))

meltA <- melt(t(A))
colnames(meltA) <- c("x", "y", "value")

ggplot(meltA, aes(x, y, fill = value) ) + geom_tile() + scale_fill_gradient(low = "darkblue", high = "green3")
ggplot(meltA, aes(x, y, fill = value) ) + geom_tile() + scale_fill_gradient(low = "darkblue", high = "greenyellow")

p <- ggplot(meltA, aes(x, y, fill = value))#
p + geom_tile() + scale_fill_gradient(low = "firebrick2", high = "yellow", 
                                      guide = guide_legend(title = "Coefficient", reverse = T),
                                      breaks = round(seq(0, max(meltA$value, na.rm = T), length.out = 10), 3)) + 
xlab("Iteration") + ylab("Feature") + theme_bw() + theme(axis.text.x=element_text(size=14),
                                                         axis.text.y=element_text(size=14), 
                                                         axis.title =element_text(size=15))


ggplot(meltA, aes(x, y, fill = value) ) + geom_tile() + scale_fill_gradient(low = "red3", high = "yellow")

ggplot(meltA, aes(x, y, fill = value) ) + geom_tile(colour = heat.colors) 

ggplot(meltA, aes(x, y, fill = value) ) + geom_tile(colour = heat.colors) 


save(A, file = "~/AExample.RData")

heatmap(t(A), Rowv = NA, Colv = NA)

load("AExample.RData")
print(levelplot(t(A), 
                xlab = list(label = "Iteration", cex = 1.5), 
                ylab = list(label = "Feature", cex = 1.5),
                col.regions = heat.colors, colorkey = list(labels = list(cex=1.5)),
                scales = list(cex = 1.5),
                panel = function(...) {
                    panel.fill(col = "black")
                    panel.levelplot(...)
                } ) )



debug(ggPlotExpression)

# flyWorm correlation

corTmp <- cor(t(flyWorm[pss.fwBC30$rowIdx, pss.fwBC30$colIdx]))
hist(corTmp[upper.tri(corTmp)], 
     prob = TRUE,
     xlab = "Correlation", 
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corFlyWormClust.pdf")

set.seed(42)
idxToSample <- setdiff(1:nrow(flyWorm), pss.fwBC30$rowIdx)
corTmp <- cor(t(flyWorm[sample(idxToSample, length(pss.fwBC30$rowIdx)), pss.fwBC30$colIdx]))
hist(corTmp[upper.tri(corTmp)], 
     prob = TRUE,
     xlab = "Correlation", 
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corFlyWormRanGene.pdf")

# perturb corr
datTmp <- sqrt(sqrt(fdr30p90))
pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.75, 0.90)
corTmp <- cor(t(datTmp[pssTmp$rowIdx, pssTmp$colIdx]))
hist(corTmp[upper.tri(corTmp)], 
     prob = TRUE,
     xlab = "Correlation", 
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corFlyPer.pdf")

set.seed(42)
idxToSample <- setdiff(1:nrow(datTmp), pssTmp$rowIdx)
corTmp <- cor(t(datTmp[sample(idxToSample, length(pss.fwBC30$rowIdx)), pssTmp$colIdx]))
hist(corTmp[upper.tri(corTmp)], 
     prob = TRUE,
     xlab = "Correlation", 
     main = "",
     xlim = c(-1, 1), col = "lightgrey")
dev.print(pdf, "../img/corFlyPerRan.pdf")

