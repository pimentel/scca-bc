# required for simulating from multivariate normal (rmvnorm)
library(mvtnorm)
library(lattice)
library(gplots)
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



covMat <- genPositiveDefMat(10, covMethod = "unifcorrmat", rangeVar = c(1,1))
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


# Progressive on genes (rows)
generateRandomBlock.progRows <- function(nrow = 4, ncol = 5, nbase = 3)
{
    ranMatrix <- sapply(1:nrow, function(rowNum)
           {
               baseNum <- rnorm(nbase + ncol, mean = 0, sd = 1)
               baseCumSum <- cumsum(baseNum)
               return(baseCumSum[(nbase+1):(nbase+ncol)])
           })
    return(list(sample = t(ranMatrix)))
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


plotParSolution <- function(parSol)
{
    par(ask = FALSE)
    A <- sapply(parSol, function (x) x$ab)
    D <- sapply(parSol, function (x) x$d)
    print(levelplot(A, xlab = "Gene", ylab = "Iteration"))
    readline("Next [Enter]\t")
    print(levelplot(D, xlab = "Condition", ylab = "Iteration"))
    readline("Next [Enter]\t")
    # heatmap(D, Colv = NA)
    heatmap(A)
    readline("Next [Enter]\t")
    heatmap(D)
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
