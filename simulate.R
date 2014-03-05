# required for simulating from multivariate normal (rmvnorm)
library(mvtnorm)
library(lattice)
library(gplots)
library(clusterGeneration) #genPositiveDefMat
# library(MASS)

##' Function to generate normal noise and Gaussian blocks
##' 
##' @param nrow integer denoting the number of rows in the matrix
##' @param ncol integer denoting the number of columns in the matrix
##' @param noiseMean a mean for the background data
##' @param noiseSD a standard deviation for the background data
##' @param clusterOptions a list of lists. See examples.
##' @return A matrix of size nrow x ncol with blocks in locations denoted by clusterOptions
generateNormal <- function(nrow = 200, ncol = 40, noiseMean = 0, noiseSD = 1, 
                           clusterOptions)
{
    # background noise
    mat <- matrix(rnorm(nrow * ncol, mean = noiseMean, sd = noiseSD), 
                  nrow = nrow, ncol = ncol)

    # TODO: make sure clusterOptions are within coords
    for (clust in clusterOptions)
    {
        clustNCols <- clust$x.end - clust$x.start + 1
        clustNRows <- clust$y.end - clust$y.start + 1
        if ("data" %in% names(clust))
        {
            mat[clust$x.start:clust$x.end, clust$y.start:clust$y.end] <- 
                clust$data
        }
    }

    return(mat)
}

# make a matrix for genes being correlated
generateRandomBlock <- function(nGenes = 4)
{
    correlations <- round(runif(round((nGenes^2-nGenes)/2), 
                                min = 0.5, max = 0.8), digits = 2)
    # correlations <- sample(c(1, -1), length(correlations), replace = TRUE) * correlations
    sig <- matrix(nrow = nGenes, ncol = nGenes)
    diag(sig) <- 1 # common variance of 1
    sig[upper.tri(sig)] <- correlations
    sig[lower.tri(sig)] <- t(sig)[lower.tri(sig)]

    # ensure it is positive definite
    # sig <- sig %*% t(sig)
    return(list(sample = rmvnorm(1, mean = rep(0, nGenes), sigma = sig, method = "svd"), sigma = sig))
}

# debug(generateRandomBlock)
# generateRandomBlock(40)

genRanPosDefMat <- function(nGenes = 10, min = 0.5, max = 0.8)
{
    correlations <- round(runif(round((nGenes^2-nGenes)/2), 
                                min = min, max = max), digits = 2)
    # correlations <- sample(c(1, -1), length(correlations), replace = TRUE) * correlations
    sig <- matrix(nrow = nGenes, ncol = nGenes)
    diag(sig) <- 1 # common variance of 1
    sig[upper.tri(sig)] <- correlations
    sig[lower.tri(sig)] <- t(sig)[lower.tri(sig)]

    # using Matrix package
    sigNear <- nearPD(sig)
    cov2cor(as.matrix(sigNear$mat))
}

mat <- genRanPosDefMat(30)
mat2 <- cor2cov(mat)


# make a matrix for genes being correlated
generateRandomBlock.posDefCov <- function(nGenes = 4)
{
    covMat <- genPositiveDefMat(nGenes, covMethod = "unifcorrmat", rangeVar = c(1,1))
    # ensure it is positive definite
    # sig <- sig %*% t(sig)
    #return(list(sample = rmvnorm(1, mean = rep(0, nGenes), sigma = covMat$Sigma, method = "svd"), sigma = sig))
    return(list(sample = mvrnorm(1, mu = rep(0, nGenes), Sigma = covMat$Sigma), sigma = covMat$Sigma))
}

ranBlock <- generateRandomBlock.posDefCov(20 * 10)
ranBlockMat <- matrix(ranBlock$sample[sample(1:length(ranBlock$sample))], nrow = 20, ncol = 10)


generateRandomBlock.progressive2 <- function(nrow = 5, ncol = 4, nbase = 15)
{
    baseSamp <- rnorm(nbase, mean = 0, sd = 1)

}

genRandomBlock.mvn <- function(nrow = 30, ncol = 20, min = 0.5, max = 0.8)
{
    covMat <- genRanPosDefMat(nrow, min, max)
    ranSamples <- mvrnorm(ncol, mu = rep(0, nrow), Sigma = covMat)
    t(ranSamples)
}

debug(genRandomBlock.mvn)
undebug(genRandomBlock.mvn)

mvnBlock <- genRandomBlock.mvn(300, 50)
mvnSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                         list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20, 
                                   data = mvnBlock)))
ps.mvnSim <- permuteSim(mvnSim)
bc.oneMvn.200.20 <- bcMultipleClusters(ps.mvnSim$permutedMat, 20, 1, 200)
bc.oneMvn.200.20.reorder <- reorderBCSol(bc.oneMvn.200.20, ps.mvnSim)

truth.oneAdd <- list(list(rowIdx = 1:25, colIdx = 1:20))
amrs.hp(truth.oneAdd, bc.oneMvn.200.20.reorder$clusters)

save(bc.oneMvn.200.20.reorder, mvnSim, file = "oneMvn.RData")

load("oneMvn.RData")
plotParSolution(bc.oneMvn.200.20.reorder$allBcSol[[1]])
A <- sapply(bc.oneMvn.200.20.reorder$allBcSol[[1]], function(x) x$ab)
D <- sapply(bc.oneMvn.200.20.reorder$allBcSol[[1]], function(x) x$d)
levelplot(abs(A), xlab = "Gene", ylab = "Condition", main = "Heatmap of AB")
levelplot(D, xlab = "Gene", ylab = "Condition", main = "Heatmap of D")


amrs.hp(truth.oneAdd, bc.oneMvn.200.15.reorder$clusters)



covMat <- genPositiveDefMat(10, covMethod = "eigen", lambdaLow = 10)
cov2cor(covMat$Sigma)

ranNorm <- mvrnorm(50, mu = rep(0, 10), Sigma = covMat$Sigma)
cor(ranNorm)

ranNorm2 <- rmvnorm(50, sigma = covMat$Sigma)
cor(ranNorm2)
cor(t(ranNorm2))
corAmongstSamples <- cor(t(ranNorm2))[upper.tri(cor(t(ranNorm2)))]
qqnorm(corAmongstSamples)
qqline(corAmongstSamples)
hist(corAmongstSamples)


ranBlock2 <- generateRandomBlock(20 * 10)
ranBlock2Mat <- matrix(ranBlock2$sample[sample(1:length(ranBlock2$sample))], nrow = 20, ncol = 10)

par(mfrow = c(2,2))
hist(abs(cor(ranBlockMat)[upper.tri(cor(ranBlockMat))]))
hist(abs(cor(t(ranBlockMat))[upper.tri(cor(t(ranBlockMat)))]))
hist(abs(cor(ranBlock2Mat)[upper.tri(cor(ranBlock2Mat))]))
hist(abs(cor(t(ranBlock2Mat))[upper.tri(cor(t(ranBlock2Mat)))]))


generateRandomBlock.posDefCov.norm <- function(nGenes = 4)
{
    covMat <- genPositiveDefMat(nGenes, covMethod = "unifcorrmat", rangeVarc(1,1))

    return(list(sample = mvrnorm(1, mu = rep(0, nGenes), Sigma = covMat$Sigma), sigma = covMat$Sigma))
}



generateRandomBlock.sum <- function(nGenes = 4, nCoefficients = 5)
{
    y <- rnorm(nCoefficients, mean = 1, sd = 1)
    # y <- rnorm(nCoefficients, mean = 0, sd = 1)
    samples <- sapply(1:nGenes, function(tmp)
                      {
                          as.numeric(y %*% runif(nCoefficients, min = 0.5, max = 1))
                      })
    list(sample = samples)
}

sumTest <- generateRandomBlock.sum(25 * 10, 5)
sumTestMat <- matrix(sumTest$sample, nrow = 25, ncol = 10)
cor(sumTestMat)


genBlock.iid.N <- function(nrow, ncol, mean, sd)
{
    matrix(rnorm(nrow * ncol, mean = mean, sd = sd), nrow = nrow, ncol = ncol)
}

genBlock.iid.N(10, 10, 4, 4)

mean(genBlock.iid.N(25, 15, 0, 1))

iidN.block <- generateNormal(nrow = 2000, ncol = 40, clusterOptions = 
               list(list(x.start = 1, x.end = 40, y.start = 1, y.end = 15, 
                         data = genBlock.iid.N(40, 15, 4, 2))))

ps.iidN.block <- permuteSim(iidN.block)
time2000 <- system.time(bc.iidN <- bcMultipleClusters(ps.iidN.block$permutedMat, 15, 1, 50))
bc.iidN.reorder <- reorderBCSol(bc.iidN, ps.iidN.block)

save(bc.iidN.reorder, file = "bc.iidN.reorder.RData")

truth.iidN <- list(list(rowIdx = 1:25, colIdx = 1:15))
amrs.hp(truth.iidN, bc.iidN.reorder$clusters)

load("bc.iidN.reorder.RData")
plotParSolution(bc.iidN.reorder$allBcSol[[1]])
A <- abs(sapply(bc.iidN.reorder$allBcSol[[1]], function (x) x$ab))
hc.A <- hclust(dist(abs(A)))
plot(hc.A)

plotParSolution(bc.iidN.reorder$allBcSol[[1]])


iidN.2block <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
               list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15, 
                         data = genBlock.iid.N(25, 15, 4, 1)),
                    list(x.start = 101, x.end = 125, y.start = 26, y.end = 40, 
                         data = genBlock.iid.N(25, 15, -4, 1))))
ps.iidN.2block <- permuteSim(iidN.2block)
bc.iidN.2block <- bcMultipleClusters(ps.iidN.2block$permutedMat, 15, 1, 50)
bc.iidN.2block.reorder <- reorderBCSol(bc.iidN.2block, ps.iidN.2block)
save(bc.iidN.2block.reorder, file = "bc.iidN.2block.reorder.RData")

load("bc.iidN.2block.reorder.RData")
plotParSolution(bc.iidN.2block.reorder$allBcSol[[1]])



amrs.hp(truth.oneAdd, bc.iidN.reorder$clusters)

levelplot(iidN.block)



# this is a basic additive model
generateRandomBlock.clust <- function(nrow = 4, ncol = 5, nCoefficients = 5, 
                                      rowMean = 2, colMean = 1)
{
    rowCenter <- rnorm(nrow, mean = 2)
    colCenter <- rnorm(ncol, mean = 1)
    samples <- matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
    samples <- t(apply(samples, 1, function(row) row + colCenter))
    samples <- apply(samples, 2, function(col) col + rowCenter)
    list(sample = samples)
}


# XXX: In this implementation, the variance explodes
generateRandomBlock.progressive <- function(nrow = 4, ncol = 5, nbase = 3)
{
    baseNum <- rnorm(nbase + nrow * ncol, mean = 0, sd = 1)
    baseCumSum <- cumsum(baseNum)
    index <- sample(nrow * ncol)
    ranMatrix <- matrix(baseCumSum[index], nrow = nrow, ncol = ncol)
    return(ranMatrix)
}
genes <- generateRandomBlock.progressive(nrow = 20, ncol = 10)

genAdditive <- function(nrow = 4, ncol = 5)
{
    rowBase <- sample(0:5, nrow, replace = TRUE)
    colBase <- sample(0:5, ncol, replace = TRUE)
    df <- matrix(rep(colBase, nrow), byrow = TRUE, nrow = nrow)
    df + rowBase + rnorm(nrow * ncol)
}

# Progressive on genes (rows)
generateRandomBlock.progRows <- function(nrow = 4, ncol = 5, nbase = 3)
{
    ranMatrix <- sapply(1:nrow, function(rowNum)
           {
               baseNum <- rnorm(nbase + ncol, mean = 0, sd = 1)
               baseCumSum <- cumsum(baseNum)
               return(baseCumSum[(nbase+1):(nbase+ncol)])
           })
    return(t(ranMatrix))
}


# TODO: simulation where one block is progressive and another is row/col additive
regBlock <- generateRandomBlock.clust(25, 15)
progRowsSamp <- generateRandomBlock.progRows(25, 15)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15,
                                       data = regBlock$sample),
                                  list(x.start = 100, x.end = 124, y.start = 26, y.end = 40,
                                       data = progRowsSamp$sample)))
bc.oneProgOneAdd.200.15 <- biclusteringPar(geneMatSim, 200, 15)
save(geneMatSim, bc.oneProgOneAdd.200.15, file = "bc.oneProgOneAdd.200.15.RData")

load("bc.oneProgOneAdd.200.15.RData")
plotParSolution(bc.oneProgOneAdd.200.15)

bc.oneProgOneAdd.200.30 <- biclusteringPar(geneMatSim, 200, 30)
save(geneMatSim, bc.oneProgOneAdd.200.30, file = "bc.oneProgOneAdd.200.30.RData")

load("bc.oneProgOneAdd.200.30.RData")
plotParSolution(bc.oneProgOneAdd.200.30)

bc.oneProgOneAdd.200.40 <- biclusteringPar(geneMatSim, 200, 40)
save(geneMatSim, bc.oneProgOneAdd.200.40, file = "bc.oneProgOneAdd.200.40.RData")

load("bc.oneProgOneAdd.200.40.RData")
plotParSolution(bc.oneProgOneAdd.200.40)

progRowsSamp <- generateRandomBlock.progRows(25, 20)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20,
                                       data = progRowsSamp$sample)))
bc.oneProg.200 <- biclusteringPar(geneMatSim, 200, 20)
save(progRowsSamp, geneMatSim, bc.oneProg.200, file = "bc.oneProg.200.RData")

plotParSolution(bc.oneProg.200)





nr <- 20
nc <- 10
mat <- matrix(runif(20 * 10, min = -1, max = 1), nrow = nr, ncol = nc)
diag(mat) <- 1



debug(generateNormal)
undebug(generateNormal)

testMat <- generateNormal(clusterOptions = 
                          list(list(x.start = 3, x.end = 10, y.start = 5, y.end = 20)))

testMat <- generateNormal(nrow = 40, ncol = 14, 
                          clusterOptions = 
                          list(list(x.start = 3, x.end = 10, y.start = 2, y.end = 3),
                               list(x.start = 13, x.end = 17, y.start = 1, y.end = 1)))
levelplot(testMat, col.regions=redgreen(75))




block <- generateRandomBlock(25)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 3, x.end = 7, y.start = 3, y.end = 7,
                                       data = block$sample)))



block <- generateRandomBlock(25*20)

block <- generateRandomBlock.sum(25*20)

block <- generateRandomBlock.clust(nrow = 25, ncol = 20)
block <- generateRandomBlock.clust(nrow = 25, ncol = 20)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20,
                                       data = block$sample)))

levelplot(geneMatSim, col.regions=redgreen(75), xlab = "Genes", ylab = "Expression")

dev.print(pdf, "../img/heatmapSim.pdf")
levelplot(scale(geneMatSim), col.regions=redgreen(75))


matplot(t(geneMatSim[1:7,]), type = "l")
dev.print(pdf, "../img/geneExpressionSim.pdf")


corCluster <- cor(t(geneMatSim[1:25, 1:20]))
corRandom <- cor(geneMatSim[1:25, 21:40])

par(mfrow = c(1,2))
hist(corCluster[upper.tri(corCluster)])
hist(corRandom[upper.tri(corRandom)])

dev.print(pdf, "../img/additiveModel.pdf")



measureCorr <- function(data, bcSol)
{
    sapply(1:length(bcSol$d), function(it)
           {
               dataX <- as.matrix(data[bcSol$splitIdx[[it]][[1]],])
               dataY <- as.matrix(data[bcSol$splitIdx[[it]][[2]],])
               lhs <- bcSol$a[[it]] %*% dataX %*% diag(bcSol$d[[it]])
               rhs <- bcSol$b[[it]] %*% dataY %*% diag(bcSol$d[[it]])
               lhs %*% t(rhs)
           })
}


measureCorrPar <- function(data, bcSol)
{
    sapply(bcSol, function(it)
           {
               dataX <- as.matrix(data[it$splitIdx[[1]],])
               dataY <- as.matrix(data[it$splitIdx[[2]],])
               lhs <- it$a %*% dataX %*% diag(it$d)
               rhs <- it$b %*% dataY %*% diag(it$d)
               lhs %*% t(rhs)
           })
}

debug(measureCorr)
measureCorr(scale(geneMatSim), bcSimSol)

bcSimSol <- biclustering(scale(geneMatSim), 50)

bcSimSol <- biclustering((geneMatSim), 50)

bcSimSolPar <- biclusteringPar(scale(geneMatSim), 400)
bcSimSolPar.noScale <- biclusteringPar(geneMatSim, 500)
bcSimSolPar.noScale.800 <- biclusteringPar(geneMatSim, 800)
bcSimSolPar.noScale.1000 <- biclusteringPar(geneMatSim, 1000)

bcSimSolPar.noScale.200 <- biclusteringPar(geneMatSim, 200)

block <- generateRandomBlock.clust(nrow = 25, ncol = 20)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20,
                                       data = block$sample)))
bcSimSolPar.200 <- biclusteringPar(geneMatSim, 200)
save(block, geneMatSim, bcSimSolPar.200, file = "bcSimSolPar.200.RData")

# no restriction on lam
block <- generateRandomBlock.clust(nrow = 25, ncol = 20)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20,
                                       data = block$sample)))
bcSimSolPar.noLam.200 <- biclusteringPar(geneMatSim, 200, 40)
save(block, geneMatSim, bcSimSolPar.noLam.200, file = "bcSimSolPar.noLam.200.RData")

# this is one smaller cluster (less conditions)
block <- generateRandomBlock.clust(nrow = 20, ncol = 10)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 20, y.start = 1, y.end = 10,
                                       data = block$sample)))
bcSimSolPar.oneSmallCluster.200 <- biclusteringPar(geneMatSim, 200)
save(block, geneMatSim, bcSimSolPar.oneSmallCluster.200, file = "bcSimSolPar.oneSmallCluster.200.RData")

load("bcSimSolPar.oneSmallCluster.200.RData")

# Testing two clusters of size 25 x 20... different means
block1 <- generateRandomBlock.clust(nrow = 25, ncol = 20)
block2 <- generateRandomBlock.clust(nrow = 25, ncol = 20, rowMean = 3, colMean = 3)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 20,
                                       data = block1$sample),
                                  list(x.start = 100, x.end = 124, y.start = 21, y.end = 40,
                                       data = block2$sample)))
bcSimSolPar.twoSmall.200 <- biclusteringPar(geneMatSim, 200)
save(block1, block2, geneMatSim, bcSimSolPar.twoSmall.200, file = "bcSimSolPar.twoSmall.200.RData")

# note: changed the value of lambda so that all conditions could be used
bcSimSolPar.twoSmall.DblLam.200 <- biclusteringPar(geneMatSim, 200)
save(block1, block2, geneMatSim, bcSimSolPar.twoSmall.DblLam.200, file = "bcSimSolPar.twoSmall.DblLam.200.RData")


plotParSolution <- function(parSol, fBase = NA, rowNames = NA,
                            colNames = NA)
{
    par(ask = FALSE)
    A <- abs(sapply(parSol, function (x) x$ab))
    D <- sapply(parSol, function (x) x$d)
    dRot <- 0
    fSize <- 1.5
    if (!is.na(rowNames)[1])
    {
        rownames(A) <- rowNames
    }
    if (!is.na(colNames)[1])
    {
        rownames(D) <- colNames
        dRot <- 45
        fSize = 1
    }
    # print(levelplot(t(A), xlab = "Gene", ylab = "Iteration",
    print(levelplot(t(A), 
                    xlab = list(label = "Iteration", cex = 1.5), 
                    ylab = list(label = "Feature", cex = 1.5),
                    col.regions = heat.colors, colorkey = list(labels = list(cex=1.5)),
                    scales = list(cex = 1.5),
                    panel = function(...) {
                        panel.fill(col = "black")
                        panel.levelplot(...)
                    } ) )
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        dev.print(pdf, paste("../img/", fBase, "Feature.pdf", sep = ""))
    print(levelplot(D, 
                    xlab = list(label = "Condition", cex = 1.5), 
                    ylab = list(label = "Iteration", cex = 1.5),
                    col.regions = heat.colors, colorkey = list(labels = list(cex = 1.5)),
                    scales = list(cex = 1.5, x = list(cex = fSize, rot = dRot))
                    ))
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        dev.print(pdf, paste("../img/", fBase, "Cond.pdf", sep = ""))
    # readline("Next [Enter]\t")
    # heatmap(A, Colv = NA)
    # readline("Next [Enter]\t")
    # heatmap(D, Colv = NA)
    return(NULL)
}

plotParSolution(bcSimSolPar.oneSmallCluster.200)
plotParSolution(bcSimSolPar.noScale.200)
plotParSolution(bcSimSolPar.twoSmall.200)
plotParSolution(bcSimSolPar.twoSmall.DblLam.200)
plotParSolution(bcSimSolPar.noLam.200)

levelplot(geneMatSim)

bcSimSolPar.scaleInternal.1000 <- biclusteringPar(geneMatSim, 1000)

save(block, geneMatSim, bcSimSolPar, file = "bcSimSolPar.RData")
save(bcSimSolPar.noScale, file = "bcSimSolPar.noScale.RData")
save(bcSimSolPar.noScale.1000, file = "bcSimSolPar.noScale.1000.RData")
save(block, geneMatSim, bcSimSolPar.noScale.200, file = "bcSimSolPar.noScale.200.RData")



# load("bcSimSolPar.RData")
# load("bcSimSolPar.noScale.RData")
# load("bcSimSolPar.noScale.800.RData")
# load("bcSimSolPar.noScale.200.RData")

load("bcSimSolPar.200.RData", verbose = TRUE)
load("bcSimSolPar.noLam.200.RData", verbose = TRUE)
load("bcSimSolPar.twoSmall.200.RData", verbose = TRUE)
load("bcSimSolPar.twoSmall.DblLam.200.RData")
load("bc.oneProg.200.RData")

measureCorrPar(scale(geneMatSim), bcSimSol)

A <- sapply(bcSimSolPar, function(x) x$ab)
D <- sapply(bcSimSolPar, function(x) x$d)

A <- sapply(bcSimSolPar.noScale, function(x) x$ab)
D <- sapply(bcSimSolPar.noScale, function(x) x$d)

A <- sapply(bcSimSolPar.noScale.800, function(x) x$ab)
D <- sapply(bcSimSolPar.noScale.800, function(x) x$d)

A <- sapply(bcSimSolPar.noScale.200, function(x) x$ab)
D <- sapply(bcSimSolPar.noScale.200, function(x) x$d)

levelplot(A)
levelplot(D)

par(mfrow = c(1,2))
hist(D[1:20,])
hist(D[21:40,])

par(mfrow = c(1,2))
hist(A[1:25,])
hist(A[25:nrow(A),])

apply(D, 1, mean)
apply(A, 1, mean)

heatmap(A)

heatmap(D)


sapply(1:length(bcSimSol$a), function(it) {
       bcSimSol$ab[[it]] %*% geneMatSim %*% bcSimSol$d[[it]]
                             })



# once you have a solution, attempt to extract a bicluster from it
A <- sapply(bcSimSolPar.200, function (x) x$ab)
D <- sapply(bcSimSolPar.200, function (x) x$d)

A <- sapply(bc.oneProgOneAdd.200.15, function (x) x$ab)
D <- sapply(bc.oneProgOneAdd.200.15, function (x) x$d)

absA <- abs(A)
distA <- as.matrix(dist(A))
kmeans(absA, 2)
kmeans(D, 2)

aHclust <- hclust(dist(absA))
cutree(aHclust, 2)


load("~/tmp/bcSession.RData")


# writing function to permute rows and columns before sending it off to BC
permuteSim <- function(simDf)
{
    rowOrder <- sample(nrow(simDf))
    colOrder <- sample(ncol(simDf))

    return(list(rowOrder = rowOrder, colOrder = colOrder, permutedMat
                = simDf[rowOrder, colOrder], mat = simDf))
}

# test using permuted sim
addBlock <- generateRandomBlock.clust(25, 15)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15,
                                       data = addBlock$sample)))
ps.geneMatSim <- permuteSim(geneMatSim)
bc.oneAdd.200.15 <- bcMultipleClusters(ps.geneMatSim$permutedMat, 15, 1, 200)

truth.oneAdd <- list(list(rowIdx = 1:25, colIdx = 1:15))

bc.oneAdd.200.15.reorder <- reorderBCSol(bc.oneAdd.200.15, ps.geneMatSim)

amrs(truth.oneAdd, bc.oneAdd.200.15.reorder$clusters)
amrs.hp(truth.oneAdd, bc.oneAdd.200.15.reorder$clusters)
save(bc.oneAdd.200.15, bc.oneAdd.200.15.reorder, ps.geneMatSim, geneMatSim, 
     file = "permuteSimTest.RData")

load("permuteSimTest.RData")

plotParSolution(bc.oneAdd.200.15$allBcSol[[1]])
plotParSolution(bc.oneAdd.200.15.reorder$allBcSol[[1]])
A <- abs(sapply(bc.oneAdd.200.15.reorder$allBcSol[[1]], function(x) x$ab))
cutree(hclust(dist(abs(A))), 2)

cutHCluster(A)


# check performance of hierarchical clustering on additive model
hcTestAdd <- lapply(1:50, function(it)
                    {
                        addBlock <- generateRandomBlock.clust(25, 15)
                        geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                                                     list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15,
                                                               data = addBlock$sample)))
                        ps.geneMatSim <- permuteSim(geneMatSim)
                        bc <- bcMultipleClusters(ps.geneMatSim$permutedMat, 15, 1, 100)

                        bc.reorder <- reorderBCSol(bc, ps.geneMatSim)
                        return(bc.reorder)
                    })

sapply(hcTestAdd, function(it) amrs(truth.oneAdd, it$clusters))
mean(sapply(hcTestAdd, function(it) amrs.hp(truth.oneAdd, it$clusters)))

# check the performance of kmeans clustering on additive model

# check performance of hierarchical clustering on additive model
kTestAdd <- lapply(1:50, function(it)
                    {
                        addBlock <- generateRandomBlock.clust(25, 15)
                        geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                                                     list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15,
                                                               data = addBlock$sample)))
                        ps.geneMatSim <- permuteSim(geneMatSim)
                        bc <- bcMultipleClusters.kmeans(ps.geneMatSim$permutedMat, 15, 1, 100)

                        bc.reorder <- reorderBCSol(bc, ps.geneMatSim)
                        return(bc.reorder)
                    })

mean(sapply(kTestAdd, function(it) amrs(truth.oneAdd, it$clusters)))
mean(sapply(kTestAdd, function(it) amrs.hp(truth.oneAdd, it$clusters)))

reorderBCSol <- function(bcSol, ps)
{
    for (bcIt in 1:length(bcSol$allBcSol))
    {
        bcSol$allBcSol[[bcIt]] <- lapply(bcSol$allBcSol[[bcIt]], function(aSol) 
                                         {
                                             aSol$ab <- aSol$ab[order(ps$rowOrder)]
                                             aSol$d <- aSol$d[order(ps$colOrder)]
                                             return(aSol)
                                         })
    }

    bcSol$clusters <- lapply(bcSol$clusters, function(clust)
                             {
                                 clust$rowIdx <- ps$rowOrder[clust$rowIdx]
                                 clust$colIdx <- ps$colOrder[clust$colIdx]
                                 return(clust)
                             })

    return(bcSol)
}

reorganizeSim <- function(psRes)
{
    psRes$permutedMat[order(psRes$rowOrder), order(psRes$colOrder)]
}

# testing permuteSim
mat <- matrix(1:20, ncol = 5, byrow = TRUE)
psMat <- permuteSim(mat)
reorganizeSim(psMat)

par(mfrow = c(2, 1))
levelplot(geneMatSim)
levelplot(permuteSim(geneMatSim)$permutedMat)


block1 <- generateRandomBlock.clust()



# testing finding multiple clusters
regBlock <- generateRandomBlock.clust(25, 15)
progRowsSamp <- generateRandomBlock.progRows(25, 15)
geneMatSim <- generateNormal(nrow = 220, ncol = 40, clusterOptions = 
                             list(list(x.start = 1, x.end = 25, y.start = 1, y.end = 15,
                                       data = regBlock$sample),
                                  list(x.start = 100, x.end = 124, y.start = 26, y.end = 40,
                                       data = progRowsSamp$sample)))
ps.geneMatSim <- permuteSim(geneMatSim)
bc.oneProgOneAdd.200.15 <- bcMultipleClusters(ps.geneMatSim$permutedMat, 15, 2, 200)

bc.oneProgOneAdd.200.15.reorder <- reorderBCSol(bc.oneProgOneAdd.200.15, ps.geneMatSim) 
save(geneMatSim, ps.geneMatSim, bc.oneProgOneAdd.200.15, bc.oneProgOneAdd.200.15.reorder,
     file = "bc.oneProgOneAdd.200.15.RData")

load("bc.oneProgOneAdd.200.15.RData")
plotParSolution(bc.oneProgOneAdd.200.15.reorder$allBcSol[[2]])

bc.oneProgOneAdd.200.15.reorder$allBcSol[[2]]$clusters


# evaluate how good the simulation is doing
# truth and pred are lists of lists. Each index in the list corresponds to a list with:
# G: the set of genes
# C: the set of conditions
amrs <- function(truthList, predList)
{
    mean(sapply(truthList, function(truth)
           {
               max(sapply(predList, function(prediction)
                      {
                          numerator <- length(intersect(truth$rowIdx, prediction$rowIdx)) 
                          denominator <- length(union(truth$rowIdx, prediction$rowIdx))
                          numerator / denominator
                      }))
           }))
}


# check how well you do on every cell
amrs.hp <- function(truthList, predList)
{
    expandCells <- function(aList)
    {
        aCat <- expand.grid(aList$rowIdx, aList$colIdx)
        paste(aCat[,1], aCat[,2], sep = ".")
    }
    mean(sapply(truthList, function(truth)
           {
               truthCat <- expandCells(truth)
               max(sapply(predList, function(prediction)
                      {
                          predCat <- expandCells(prediction)
                          numerator <- length(intersect(truthCat, predCat)) 
                          denominator <- length(union(truthCat, predCat))
                          numerator / denominator
                      }))
           }))
}


# check how well you do on every cell
amrs.eachClust <- function(truthList, predList)
{
    expandCells <- function(aList)
    {
        aCat <- expand.grid(aList$rowIdx, aList$colIdx)
        paste(aCat[,1], aCat[,2], sep = ".")
    }
    (sapply(truthList, function(truth)
           {
               truthCat <- expandCells(truth)
               max(sapply(predList, function(prediction)
                      {
                          predCat <- expandCells(prediction)
                          numerator <- length(intersect(truthCat, predCat)) 
                          denominator <- length(union(truthCat, predCat))
                          numerator / denominator
                      }))
           }))
}


# check how well you do on genes then on conditions
amrs.hp2 <- function(truthList, predList)
{
    mean(sapply(truthList, function(truth)
           {
               max(sapply(predList, function(prediction)
                      {
                          #predCat <- expandCells(prediction)
                          numerator <- length(intersect(truth$rowIdx, prediction$rowIdx))  
                          numerator <- numerator + length(intersect(truth$colIdx, prediction$colIdx))
                          denominator <- length(union(truth$rowIdx, prediction$rowIdx))
                          denominator <- denominator + length(union(truth$colIdx, prediction$colIdx))
                          numerator / denominator
                      }))
           }))
}

t.amrs.truth <- list(list(rowIdx = 1:10, colIdx = 2:4),
                     list(rowIdx = 5:20, colIdx = 3:10))
t.amrs.pred <- list(list(rowIdx = 1:9, colIdx = 2:4),
                    list(rowIdx = c(5, 7, 14:18), colIdx = 3:9))

amrs(t.amrs.truth, t.amrs.pred)
amrs(t.amrs.truth, t.amrs.truth)
amrs(t.amrs.pred, t.amrs.truth)


amrs.hp(t.amrs.truth, t.amrs.pred)
amrs.hp(t.amrs.truth, t.amrs.truth)
amrs.hp(t.amrs.pred, t.amrs.truth)

amrs.hp2(t.amrs.truth, t.amrs.pred)
amrs.hp2(t.amrs.truth, t.amrs.truth)
amrs.hp2(t.amrs.pred, t.amrs.truth)



# testing sub sampling
iidN.block <- generateNormal(nrow = 300, ncol = 40, clusterOptions = 
               list(list(x.start = 1, x.end = 40, y.start = 1, y.end = 15, 
                         data = genBlock.iid.N(40, 15, 4, 2))))

levelplot(iidN.block[sort(sample(nrow(iidN.block), round(nrow(iidN.block)*0.6) )),])
iid.subSample <- bcSubSamplePar(iidN.block, 100, 15, 0.6)
save(iid.subSample, file = "iid.subSample.6.RData")

load("iid.subSample.6.RData")
plotParSolution(iid.subSample)

debug(postSubSample)
post.iid.subSample <- postSubSample(iid.subSample, 0.95)
post.iid.subSample <- postSubSample.percent(iid.subSample, 0.95, 0.7)
apply(post.iid.subSample$ab, 2, mean)
plot(post.iid.subSample$ab)
levelplot(post.iid.subSample$ab)
levelplot(post.iid.subSample$d)

kClust <- kmeans(post.iid.subSample$ab, 3)

cutree(hclust(dist(post.iid.subSample$ab)), 3)

hPclust <- cutree(hclust(dist(post.iid.subSample$ab)), 3)
hPclust

truth.iid <- list(list(rowIdx = 1:40, colIdx = 1:15))
amrs.hp(truth.iid, post.kmeans(post.iid.subSample))
amrs.hp(truth.iid, post.hclust(post.iid.subSample))

post.kmeans(post.iid.subSample)



iidN.block <- generateNormal(nrow = 180, ncol = 40, clusterOptions = 
               list(list(x.start = 1, x.end = 24, y.start = 1, y.end = 15, 
                         data = genBlock.iid.N(24, 15, 4, 2))))

iid.fullSample <- biclusteringPar(iidN.block, 100, 15)
save(iid.fullSample, file = "iid.fullSample.RData")
load("iid.fullSample.RData")
plotParSolution(iid.fullSample)

iid.subSample <- bcSubSamplePar(iidN.block, 15, 2, 0.6)


debugonce(ggPlotParSolution)

ggPlotParSolution <- function(parSol, fBase = NA, rowNames = NA,
                              colNames = NA)
{
    par(ask = FALSE)
    A <- abs(sapply(parSol, function (x) x$ab))
    D <- sapply(parSol, function (x) x$d)
    dRot <- 0
    dxSize <- 14
    fSize <- 1.5
    if (!is.na(rowNames)[1])
    {
        rownames(A) <- rowNames
    }
    if (!is.na(colNames)[1])
    {
        rownames(D) <- colNames
        dRot <- 90
        dxSize <- 10
        fSize = 1
    }
    meltA <- melt(t(A))
    colnames(meltA) <- c("x", "y", "value")
    breaksA <- round(seq(0, max(meltA$value, na.rm = T), length.out = 10), 3)
    p <- ggplot(meltA, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradient(low = "firebrick2", high = "yellow", 
                                               guide = guide_legend(title = "Coefficient", reverse = T), 
                                               breaks = breaksA) 

    p <- p + xlab("Iteration") + ylab("Feature") + theme_bw() + theme(axis.text.x=element_text(size=14),
                                                                      axis.text.y=element_text(size=14), 
                                                                      axis.title=element_text(size=15),
                                                                      legend.text=element_text(size=14),
                                                                      legend.title=element_text(size=14))
    print(p)
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        ggsave(paste("../img/gg", fBase, "Feature.pdf", sep = ""), 
               width = 6.62, height = 12.8)
    meltD <- melt(D)
    colnames(meltD) <- c("x", "y", "value")
    breaksD <- round(seq(0, max(meltD$value, na.rm = T), length.out = 10), 1)
    p <- ggplot(meltD, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradient(low = "firebrick2", high = "yellow", 
                                               guide = guide_legend(title = "Coefficient", reverse = T), 
                                               breaks = breaksD) 
    p <- p + xlab("Condition") + ylab("Iteration") + theme_bw() + theme(axis.text.x=element_text(size=dxSize,
                                                                                                 angle = dRot),
                                                    axis.text.y=element_text(size=14), 
                                                    axis.title =element_text(size=15),
                                                    legend.text=element_text(size=14),
                                                    legend.title=element_text(size=14))
    print(p)
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        ggsave(paste("../img/gg", fBase, "Cond.pdf", sep = ""),
               width = 6.62, height = 12.8)

    # Saving 6.62 x 12.8 in image
    return(NULL)
}

ggPlotParSolution(nearMedMvn, )
