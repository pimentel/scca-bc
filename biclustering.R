# Author: Harold Pimentel
# Contact: pimentel@cs.berkeley.edu

# Lee et. al. package. Provides 'scca()'
library(scca)
library(MASS)
library(Matrix)
library(multicore)

debug(biClustMax.optim)
undebug(biClustMax.optim)

debug(NIPALS.sparse)
undebug(NIPALS.sparse)

debug(opt.cv.alt)
undebug(opt.cv.alt)

debug(scca)
undebug(scca)

##' Randomly splits indices if there are an even number
##' 
##' @param n the number of indices
##' @returns a list with two disjoint sets of indices from 1 to n.
splitEvenly <- function(n)
{
    if ( n %% 2 != 0)
    {
        stop("Not an even number of genes!")
    }
    sets <- split(sample.int(n), 1:2)
    return(sets)
}

##' One iteration of biclustering maximization. Takes the current values of 
##' d, a, b and if successful returns one iteration
##' 
##' @param xDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param yDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param d vector of dimension k
##' @param a vector of dimension n
##' @param b vector of dimension n
##' @param lam maximum cluster size
##' @return A list containing maximum values of a, b, d
# uses constrOptim
biClustMax.optim <- function(xDf, yDf, d.start, a, b, lam, verbose = TRUE, optim.max = 4, lam.lwr = 3.5,
                             clustOptions = list())
{

    # iterate until we've hit maximum iterations or QP has converged
    if (verbose)
        cat("Computing tilde matrices\n")
    xTilde <- t(xDf %*% diag(d.start))
    yTilde <- t(yDf %*% diag(d.start))

    if (verbose)
        cat("Computing SCCA component\n")

    # First maximize a and b using SCCA
    pen <- "LASSO"
    if (!is.null(clustOptions$penalty))
        pen <- clustOptions$penalty

    lamx <- c(1,2,3)
    if (!is.null(clustOptions$lamx))
        lamx <- clustOptions$lamx
    cv <- "CV.alt"
    if (!is.null(clustOptions$cv))
        cv <- clustOptions$cv

    if (verbose)
        cat("Using penalty: ", pen, "\n")

    sccaRes <- scca(as.matrix(xTilde), as.matrix(yTilde), 
                    nc = 1, 
                    penalty = pen,
                    lamx = lamx,
                    lamy = lamx,
                    tuning = cv,
                    center = TRUE, scale = TRUE)
    a <- as.numeric(sccaRes$A)
    b <- as.numeric(sccaRes$B)

    # XXX: Possibly remove... Unsure if I should scale in this fashion
    xDf <- scale(xDf)
    yDf <- scale(yDf)

    # Now, maximize D using quadratic programming
    # q <- - 2 * as.numeric((a %*% xDf) * (b %*% yDf))
    q <- -as.numeric((a %*% xDf) * (b %*% yDf))

    if (verbose)
        cat("Solving optimization of d\n")

    # lower bound is zero
    constrMat<- diag(length(q))
    # upper bound is one
    constrMat <- rbind(constrMat, -diag(length(q)))
    # regularize (sum(d) <= lam)
    constrMat <- rbind(constrMat, rep.int(-1, length(q)))
    # regularize (sum(d) >= lam.lwr)
    constrMat <- rbind(constrMat, rep.int(1, length(q)))
    constrLimit <- c(rep(0, length(q)), 
                     rep(-1, length(q)),
                     -lam, lam.lwr)
    objectiveFn <- function(d, qVec)
    {
        # sum(0.5 * d^2 * qVec)
        sum(d^2 * qVec)
    }

    optimIt <- 0
    repeat {
        optimRes <- constrOptim(d.start, objectiveFn, grad = NULL,
                                control = list(maxit = as.integer(1000000000)),
                                # control = list(maxit = 10000),
                                ui = constrMat, ci = constrLimit,
                                qVec = q
                                )

        if (optimRes$convergence != 0)
        {
            warning("Error with convergence.\n", "\tError code: ", 
                    optimRes$convergence, "\tMessage: ", optimRes$message,
                    immediate. = TRUE)
            if (optimIt >= optim.max)
            {
                warning("*****Reached maximum number of iterations of constrOptim. Exiting",
                        immediate. = TRUE)
                break
            }
            optimIt <- optimIt + 1
            d.start <- optimRes$par
            cat("\t\tRestarting optimization at new point.\n")
            cat("\t\tconstrOptim iteration: ", optimIt + 1, "\n")
        }
        else
        {
            break
        }

    }

    return(list(a = a, b = b, d = optimRes$par, q, sccaLam = sccaRes$lambda))
}


# TODO: Replace the regular distance in maximizeOneSplit with pDiff and see how
# affects convergence
pDiff <- function(old, new)
{
    num <- dist(rbind(old, new))[1]
    denom <- sqrt(sum(old^2))
    num / denom
}

##' Driver function to find solution for 
maximizeOneSplit <- function(geneDf, lam, epsA = 0.1, epsB = 0.1, epsD = 0.1, maxIt = 100, lam.lwr, 
                             clustOptions = list())
{
    splitIdx <- splitEvenly(nrow(geneDf))
    xDf <- geneDf[splitIdx[[1]], ]
    yDf <- geneDf[splitIdx[[2]], ]

    # TODO: implement cross validation for lambda
    # lam <- round(min(dim(xDf)) * .3)
    # lam <- 20 * 2


    # randomly initialize d
    # FIXME: If there are lots of conditions, possible that won't be able to
    # find d that satisfy the constraint (if lam is sufficiently small. Come up
    # with different way to randomly assign values
    d <- 0
    minUnif <- 0.05
    maxUnif <- 1
    repeat {
        d <- runif(ncol(xDf), min = minUnif, max = maxUnif)
        sumD <- sum(d)
        if (sumD <= lam & sumD >= lam.lwr)
            break
        else if (sumD < lam.lwr)
            minUnif <- min(minUnif + 0.03, maxUnif - 0.03)
        else
            maxUnif <- max(minUnif + 0.04, maxUnif - 0.04)
    }
    if (sum(d) < lam.lwr)
    {
        cat("ERROR! sum(d) < lam.lwr\n")
        cat("ERROR! sum(d) < lam.lwr\n")
        cat("ERROR! sum(d) < lam.lwr\n")
        cat("ERROR! sum(d) < lam.lwr\n")
        cat("ERROR! sum(d) < lam.lwr\n")
        cat("ERROR! sum(d) < lam.lwr\n")
    }
    cat("sum(d)", sum(d), "\n") 
    # NB: values of a and b are not important currently since evaluated after d is set.
    # If ever change the order, fix this.
    a <- b <- runif(nrow(xDf), min = -1, max = 1)
    it <- 1
    curSol <- 0
    repeat {
        cat("\tOne split iteration: ", it, "\n")
        curSol <- biClustMax.optim(xDf, yDf, d, a, b, lam, lam.lwr = lam.lwr, clustOptions = clustOptions)
        cat ("\t\tdist: ",
             dist(rbind(curSol$a, a))[1], "\t",
             dist(rbind(curSol$b, b))[1], "\t",
             dist(rbind(curSol$d, d))[1], "\n")
        if (dist(rbind(curSol$a, a))[1] < epsA &
            dist(rbind(curSol$b, b))[1] < epsB &
            dist(rbind(curSol$d, d))[1] < epsD & it > 2) # potentially remove it > 2
        {
            cat("\tConverged.\n")
            break
        }

        a <- curSol$a
        b <- curSol$b
        d <- curSol$d
        it <- it + 1
        if (it > maxIt)
        {
            cat("Warning: exceed maximum number of iterations (", maxIt, ")\n")
            break
        }
    }
    # NB: If there is something funky with the distribution of a and b, 
    # check here... though this should work correctly
    # a <- rep.int(0, length(a))
    # b <- rep.int(0, length(b))
    ab <- rep.int(0, length(b*2))
    ab[splitIdx[[1]]] <- curSol$a
    ab[splitIdx[[2]]] <- curSol$b
    # a[splitIdx[[1]]] <- curSol$a
    # b[splitIdx[[2]]] <- curSol$b

    return(list(a = a, b = b, ab = ab, d = round(curSol$d, digits = 4),
                splitIdx = splitIdx, sccaLam = curSol$sccaLam))
}


#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(cores = N).
#' @param geneDf data.frame with genes on rows and columns defining conditions
#' @param nSamples integer denoting the number of permutations to perform
#' @param lam regularization parameter for conditions (maximum number of conditions allows in a bicluster)
biclusteringPar <- function(geneDf, nSamples = 100, lam)
{
    mclapply(1:nSamples, function(it) {
             cat("Biclustering iteration: ", it, "\n")
             curSol <- maximizeOneSplit(geneDf, lam)
             return(curSol)
            })
}


bcSubSamplePar <- function(geneDf, nSamples = 100, lam, propSample = 0.6, lam.lwr = 3.5,
                           clustOptions = list())
{
    if (propSample > 1 | propSample < 0.01)
        stop("Invalid range for propSample")

    nRowsSample <- round(propSample * nrow(geneDf) )
    if (nRowsSample %% 2)
        nRowsSample <- nRowsSample + 1

    cat("Sampling ", nRowsSample, " features\n")
    mclapply(1:nSamples, function(it) {
             cat("**Biclustering iteration: ", it, "\n")
             sampIdx <- sample.int(nrow(geneDf), size = nRowsSample)
             abSol <- rep.int(NA, nrow(geneDf))
             curSol <- maximizeOneSplit(geneDf[sampIdx,], lam, lam.lwr = lam.lwr, 
                                        clustOptions = clustOptions)
             abSol[sampIdx] <- curSol$ab
             d <- curSol$d
             return(list(ab = abSol, d = d, sccaLam = curSol$sccaLam))
                           })
}

postSubSample <- function(subSampleSol, percentile = 0.75)
{
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abGt <- apply(ab, 2, function(col) {
                  # median(col, na.rm = TRUE)
                  mean(col, na.rm = TRUE)
                  # quantile(col, probs = percentile, na.rm = TRUE)
            })
    dGt <- apply(d, 2, function(col) 
                  {
                      median(col, na.rm = TRUE)
                  })
    # dGt <- t(dGt)
    return(list(ab = abGt, d = dGt))
}


postSubSample.percent <- function(subSampleSol, percentile = 0.95, eps = 0.05)
{
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- sapply(subSampleSol, function (x) x$d)
    # want to implement in top 5%, 90% it is chosen
    abGt <- apply(ab, 1, function(row) {
                  # median(col, na.rm = TRUE)
                  row >= quantile(row, probs = percentile, na.rm = TRUE)
            })
    rowIdx <- which(apply(abGt, 1, mean, na.rm = T) >= eps)
    # dGt <- apply(d, 1, function(row) 
    #               {
    #                   row >= quantile(row, probs = percentile, na.rm = TRUE)
    #               })
    # dGt <- t(dGt)
    # colIdx <- which(apply(d, 1, mean, na.rm = T) >= eps)

    
    cutHCluster <- function(mat)
    {
        hc <- hclust(dist(mat))
        hcCuts <- cutree(hc, 2)

        cutIdx <- 1
        if (mean(mat[which(hcCuts == 1), ]) < 
            mean(mat[which(hcCuts == 2), ]))
        {
            cutIdx <- 2
        }
        cutIdx <- which(hcCuts == cutIdx)

        return(cutIdx)
    }

    colIdx <- cutHCluster(d)


    return(list(rowIdx = rowIdx, colIdx = colIdx))
}

debug(postSubSample.percent)

postSubSample.percent2 <- function(subSampleSol, percentile1 = 0.95, eps1 = 0.5,
                                   percentile2 = 0.95, eps2 = 0.5)
{
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))
    # want to implement in top 5%, 90% it is chosen
    abGt <- apply(ab, 1, function(row) {
                  # median(col, na.rm = TRUE)
                  row >= quantile(row, probs = percentile1, na.rm = TRUE)
            })
    rowIdx <- which(apply(abGt, 1, mean, na.rm = T) >= eps1)

    
    dGt <- apply(d, 1, function(row) {
                 row >= quantile(row, probs = percentile2, na.rm = T)
            })
    colIdx <- which(apply(dGt, 1, mean, na.rm = T) >= eps2)

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}

postSubSample.sort <- function(subSampleSol, percentile1 = 0.95, eps1 = 0.5,
                               percentile2 = 0.95, eps2 = 0.5)
{
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))
    # want to implement in top 5%, 90% it is chosen
    abGt <- apply(ab, 1, function(row) {
                  # median(col, na.rm = TRUE)
                  row >= quantile(row, probs = percentile1, na.rm = TRUE)
                               })
    rowIdx <- which(apply(abGt, 1, mean, na.rm = T) >= eps1)
    
    dSort <- t(apply(d, 2, function(col) {sort(col)}))
    cutHCluster <- function(mat)
    {
        hc <- hclust(dist(mat))
        hcCuts <- cutree(hc, 2)

        cutIdx <- 1
        if (mean(mat[which(hcCuts == 1), ]) < 
            mean(mat[which(hcCuts == 2), ]))
        {
            cutIdx <- 2
        }
        cutIdx <- which(hcCuts == cutIdx)

        return(cutIdx)
    }
    colIdx <- cutHCluster(dSort)

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}

debug(postSubSample.sort)

postSubSample.top <- function(subSampleSol, percentile1 = 0.95, eps1 = 0.5,
                               nc)
{
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- sapply(subSampleSol, function (x) x$d)
    # want to implement in top 5%, 90% it is chosen
    abGt <- apply(ab, 1, function(row) {
                  # median(col, na.rm = TRUE)
                  row >= quantile(row, probs = percentile1, na.rm = TRUE)
                               })
    rowIdx <- which(apply(abGt, 1, mean, na.rm = T) >= eps1)
    
    dMean <- apply(d, 1, mean)
    colIdx <- order(dMean, decreasing = T)[1:nc]

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}


postSubSample.median <- function(subSampleSol, abQuant = 0.5)
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abMed <- apply(ab, 2, quantile, probs = abQuant, na.rm = T)
    dMed <- apply(d, 2, median)

    return(list(ab = abMed, d = dMed))
}


postSubSample.median.kmeans <- function(subSampleSol, abQuant = 0.5, dQuant = 0.5,
                                      seMethod = "firstSEmax")
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abQuant <- apply(ab, 2, quantile, probs = abQuant, na.rm = T)
    abSd <- apply(ab, 2, sd, na.rm = T)
    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)

    clustKmeans <- function(dat)
    {
        gap <- clusGap(as.matrix(dat), kmeans, 5)
        nk <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- 2
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }
    debugonce(clustKmeans)

    rowIdx <- clustKmeans(cbind(abQuant, abSd, abQuant / abSd))
    colIdx <- clustKmeans(cbind(dQuant, dSd))

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}


# credit to stackoverflow: 
# http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode
findMode <- function(x, na.rm = T) 
{
    if (na.rm)
        x <- x[!is.na(x)]
    ux <- unique(x)
    if (length(ux) == length(x))
        print("ERROR: No mode")
    ux[which.max(tabulate(match(x, ux)))]
}

postSubSample.median.kmeans2 <- function(subSampleSol, abQuant = 0.5, dQuant = 0.5,
                                         seMethod = "firstSEmax")
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abStats <- data.frame(med = apply(ab, 2, quantile, probs = abQuant, na.rm = T),
                          upr = apply(ab, 2, quantile, probs = 0.75, na.rm = T),
                          mode = apply(round(ab, digits = 3), 2, findMode, na.rm = T),
                          sd = apply(ab, 2, sd, na.rm = T),
                          min = apply(ab, 2, min, na.rm = T))
    abStats$medRate <- with(abStats, med / sd)
    abStats$mean <- apply(ab, 2, mean, na.rm = T)
    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)

    # abMode <- apply(round(ab, digits = 3), 2, findMode, na.rm = T)

    clustKmeans <- function(dat, minK = 2)
    {
        gap <- clusGap(as.matrix(dat), kmeans, 5)
        nk <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- minK
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }
    # debugonce(clustKmeans)

    # rowIdx <- clustKmeans(abStats, 3)
    rowIdx <- clustKmeans(with(abStats, cbind(medRate, sd)), 3)
    # clustKmeans(with(abStats, cbind(med, sd)), 2)
    # clustKmeans(with(abStats, cbind(mode, medRate)), 3)
    # clustKmeans(with(abStats, cbind(mean, medRate)), 3)
    # clustKmeans(with(abStats, cbind(mean, sd, medRate)), 3)
    # clustKmeans(with(abStats, cbind(medRate, sd)), 4)

    # mat <- with(abStats, as.matrix(cbind(medRate, sd)))
    colIdx <- clustKmeans(cbind(dQuant, dSd))

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}

vScore <- function(set1, set2)
{
    length(intersect(set1, set2)) / length(union(set1, set2))
}

vScore2 <- function(...)
{
    length(intersect(...)) / length(union(...))
}

postSubSample.tight <- function(subSampleSol, abThresh = 0.6, dQuant = 0.5)
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abStats <- data.frame(med = apply(ab, 2, median, na.rm = T),
                          upr = apply(ab, 2, quantile, probs = 0.75, na.rm = T),
                          mode = apply(round(ab, digits = 3), 2, findMode, na.rm = T),
                          sd = apply(ab, 2, sd, na.rm = T),
                          min = apply(ab, 2, min, na.rm = T))
    abStats$medRate <- with(abStats, med / sd)
    abStats$mean <- apply(ab, 2, mean, na.rm = T)

    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)


    tight <- function(dat, minK = 2, maxK = 7, thresh = abThresh)
    {
        getBestClust <- function(df, theK)
        {
            kRes <- kmeans(df, theK, nstart = 20)
            kMax <- which.max(kRes$centers[,1])

            which(kRes$cluster == kMax)

        }
        clusts <- lapply(minK:maxK, function(curK)
                         {
                             getBestClust(dat, curK)
                         })
        vs <- sapply((minK+1):maxK, function(curK)
                     {
                         vScore( clusts[[curK - minK]], clusts[[curK - minK + 1]] )
                     })
        vsGt <- min(which(vs >= thresh))
        bestK <- minK + vsGt

        cat("Choosing tight cluster: ", bestK, "\n")
        getBestClust(dat, bestK)
    }

    rowIdx <- with(abStats, tight(cbind(medRate, sd)))
    # rowIdx <- with(abStats, tight(cbind(medRate)))
    # rowIdx <- with(abStats, tight(cbind(med, sd)))

    clustKmeans <- function(dat, minK = 2)
    {
        gap <- clusGap(as.matrix(dat), kmeans, 5)
        nk <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- minK
        cat("Conditions k: ", nk, "\n")
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }

    # colIdx <- clustKmeans(cbind(dQuant, dSd))
    colIdx <- clustKmeans(cbind(dQuant))

    list(rowIdx = rowIdx, colIdx = colIdx)
}


postSubSample.pca <- function(subSampleSol, abThresh = 0.6, dQuant = 0.5)
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))


    abMed <- apply(ab, 2, median, na.rm = T)
    for (col in 1:ncol(ab))
    {
        rng <- quantile(ab[,col], probs = c(0.05, 0.95), na.rm = T)
        whichValid <- which(rng[1] <= ab[,col] & ab[,col] <= rng[2])
        nNa <- sum(is.na(ab[,col]))
        ab[is.na(ab[,col]), col] <- ab[sample(whichValid, nNa, replace = T), col]
    }


    # first, 
    pcaAB <- prcomp(ab, center = F, scale. = F)
    pcaAB$rotation <- as.data.frame(data.matrix(as.data.frame(pcaAB$rotation, stringsAsFactors = F)))
    #pcaAB$rotation <- cbind(pcaAB$rotation, inBC = c(rep("bc", 300), rep("bg", 1200)))

    pcaAB$rotation[,"median"] <- apply(ab, 2, median, na.rm = T)
    pcaAB$rotation[,"var"] <- apply(ab, 2, var, na.rm = T)

    # gg <- ggplot(pcaAB$rotation, aes(x = PC1, y = PC2, colour = inBC))
    # gg <- ggplot(pcaAB$rotation, aes(x = PC1, y = median, colour = inBC))
    # gg <- ggplot(pcaAB$rotation, aes(x = PC1, y = var, colour = inBC))
    # gg <- ggplot(pcaAB$rotation, aes(x = median, y = median, colour = inBC))
    # gg + geom_point()
    # ggsave("../img/pcaMvnNoCenterNoScale.pdf")

    # plot(pcaAB$rotation[,"PC1"], pcaAB$rotation[,"PC2"])

    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)

    # test <- tight.clust(pcaAB$rotation[,c("PC1", "PC2")], 1, 3, standardize.gene = F, top.can = 1)

    tight <- function(dat, minK = 2, maxK = 7, thresh = abThresh)
    {
        # TODO: Investigate how performs w/ and w/o scaling
        dat <- scale(dat)
        getBestClust <- function(df, theK)
        {
            kRes <- kmeans(df, theK, nstart = 20)
            # kMax <- which.max(kRes$centers[,1])
            kMax <- which.min(kRes$centers[,1])

            which(kRes$cluster == kMax)

        }
        clusts <- lapply(minK:maxK, function(curK)
                         {
                             getBestClust(dat, curK)
                         })
        vs <- sapply((minK+1):maxK, function(curK)
                     {
                         vScore( clusts[[curK - minK]], clusts[[curK - minK + 1]] )
                     })
        vsGt <- min(which(vs >= thresh))
        bestK <- minK + vsGt

        cat("Choosing tight cluster: ", bestK, "\n")
        getBestClust(dat, bestK)
    }

    # tight2 <- function(dat)
    # {
    #     tc <- tight.clust(dat, 1, 3, standardize.gene = F, top.can = 1)
    # }

    # rowIdx <- with(abStats, tight(cbind(medRate, sd)))
    # rowIdx <- with(abStats, tight(cbind(medRate)))
    # rowIdx <- with(abStats, tight(cbind(med, sd)))

    # rowIdx <- with(pcaAB$rotation, tight(cbind(median, PC1, PC2)))
    # rowIdx <- with(pcaAB$rotation, tight(cbind(median, PC1)))
    rowIdx <- with(pcaAB$rotation, tight((PC1)))

    clustKmeans <- function(dat, minK = 2)
    {
        gap <- clusGap(as.matrix(dat), kmeans, 5)
        nk <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- minK
        cat("Conditions k: ", nk, "\n")
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }

    # colIdx <- clustKmeans(cbind(dQuant, dSd))
    colIdx <- clustKmeans(cbind(dQuant))

    return(list(abDat = pcaAB$rotation, cluster = list(rowIdx = rowIdx, colIdx = colIdx)))
}



postSubSample.mean.kmeans <- function(subSampleSol, abQuant = 0.5, dQuant = 0.5,
                                      seMethod = "firstSEmax")
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    abQuant <- apply(ab, 2, mean, na.rm = T)
    dQuant <- apply(d, 2, mean,na.rm = T)

    clustKmeans <- function(dat)
    {
        gap <- clusGap(as.matrix(dat), kmeans, 5)
        nk <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = seMethod)
        if (nk == 1)
            nk <- 2
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers)

        which(kRes$cluster == kMax)
    }

    rowIdx <- clustKmeans(abQuant)
    colIdx <- clustKmeans(dQuant)

    return(list(rowIdx = rowIdx, colIdx = colIdx))
}



debug(postSubSample.top)

post.hclust <- function(postSample, nAB = 3, nD = 2)
{
    clustAB <- cutree(hclust(dist(postSample$ab)), nAB)
    clustABCenter <- sapply(unique(clustAB), 
                            function(it) median(postSample$ab[which(clustAB == it)]))
    clustABMax <- which.max(clustABCenter)
    rowIdx <- which(clustAB == clustABMax)

    clustD <- cutree(hclust(dist(postSample$d)), nD)
    clustDCenter <- sapply(unique(clustD), 
                           function(it) median(postSample$d[which(clustD == it)]))
    clustDMax <- which.max(clustDCenter)
    colIdx <- which(clustD == clustDMax)
    list(list(rowIdx = rowIdx, colIdx = colIdx))
}

post.kmeans <- function(postSample, nAB = 3, nD = 2)
{
    kClustAB <- kmeans(postSample$ab, nAB)
    kClustABMax <- which.max(kClustAB$centers)
    rowIdx <- which(kClustAB$cluster == kClustABMax)

    clustD <- kmeans(postSubSample$d, nD)
    clustDMax <- which.max(clustD$center)
    colIdx <- which(clustD$cluster == clustDMax)

    list(list(rowIdx = rowIdx, colIdx = colIdx))
}


bcMultipleClusters.subSample <- function(geneDf, lam, nClusters, nSamples = 100, 
                                         propSample = 0.6)
{
    allBcSol <- list()
    clusters <- list()
    for (clust in 1:nClusters)
    {
        bcSol <- bcSubSamplePar(geneDf, nSamples, lam, propSample)
        allBcSol <- append(allBcSol, list(bcSol))

        postSol <- postSubSample.percent(bcSol, 0.9, 0.5)

        rowIdx <- postSol$rowIdx
        colIdx <- postSol$colIdx

        clusters <- append(clusters, list(list(rowIdx = rowIdx, colIdx = colIdx)))

        # given a list of rows and columns, converts pairwise combinations into
        # 1D index for matrix
        get1DIdx <- function(rows, cols) 
        {
            allPairs <- expand.grid(rows, cols)
            nrow(geneDf) * (allPairs[, 2] - 1) + allPairs[, 1]
        }

        # mask each gene with random values from the rest of the matrix also
        # considered doing this with JUST the other values of the gene... need
        # to experiment with that
        clustIdx <- get1DIdx(rowIdx, colIdx)
        geneDf[clustIdx] <- sample(geneDf[-clustIdx], length(clustIdx))
        # temporarily sample rnorm to debug FP issue
        # geneDf[clustIdx] <- rnorm(length(clustIdx))
    }

    return(list(allBcSol = allBcSol, clusters = clusters))
}




bcMultipleClusters <- function(geneDf, lam, nClusters, nSamples = 100)
{
    allBcSol <- list()
    clusters <- list()
    for (clust in 1:nClusters)
    {
        bcSol <- biclusteringPar(geneDf, nSamples, lam)
        allBcSol <- append(allBcSol, list(bcSol))

        # Take abs() of A since depending on partition, correlation could be in
        # opposite direction
        A <- abs(sapply(bcSol, function(x) x$ab))
        D <- sapply(bcSol, function(x) x$d)

        # find features/conditions which are clustered
        # to find the cluster, find cut with fewest elements...  

        # XXX: this is not necessarily correct. For example, it is *possible*
        # that a cluster might contain more than half of the genes and only
        # a few conditions.  Think about a "better" way to do this... i.e.
        # perhaps a matrix multiplication

        # returns which indices correspond to a cluster
        # TODO: fix this to get the cluster with the most non-zero entries
        cutHCluster <- function(mat)
        {
            hc <- hclust(dist(mat))
            hcCuts <- cutree(hc, 2)

            # cutIdx <- 1
            # if (mean(hcCuts == 1) > 0.5)
            # {
            #     cutIdx <- 2
            # }
            # cutIdx <- which(hcCuts == cutIdx)

            cutIdx <- 1
            if (mean(mat[which(hcCuts == 1), ]) < 
                mean(mat[which(hcCuts == 2), ]))
            {
                cutIdx <- 2
            }
            cutIdx <- which(hcCuts == cutIdx)

            return(cutIdx)
        }

        rowIdx <- cutHCluster(A)
        colIdx <- cutHCluster(D)

        clusters <- append(clusters, list(list(rowIdx = rowIdx, colIdx = colIdx)))

        # given a list of rows and columns, converts pairwise combinations into
        # 1D index for matrix
        get1DIdx <- function(rows, cols) 
        {
            allPairs <- expand.grid(rows, cols)
            nrow(geneDf) * (allPairs[, 2] - 1) + allPairs[, 1]
        }

        # mask each gene with random values from the rest of the matrix also
        # considered doing this with JUST the other values of the gene... need
        # to experiment with that
        clustIdx <- get1DIdx(rowIdx, colIdx)
        geneDf[clustIdx] <- sample(geneDf[-clustIdx], length(clustIdx))
        # temporarily sample rnorm to debug FP issue
        # geneDf[clustIdx] <- rnorm(length(clustIdx))
    }

    return(list(allBcSol = allBcSol, clusters = clusters))
}




bcMultipleClusters.kmeans <- function(geneDf, lam, nClusters, nSamples = 100)
{
    allBcSol <- list()
    clusters <- list()
    for (clust in 1:nClusters)
    {
        bcSol <- biclusteringPar(geneDf, nSamples, lam)
        allBcSol <- append(allBcSol, list(bcSol))

        # Take abs() of A since depending on partition, correlation could be in
        # opposite direction
        A <- abs(sapply(bcSol, function(x) x$ab))
        D <- sapply(bcSol, function(x) x$d)

        # find features/conditions which are clustered
        # to find the cluster, find cut with fewest elements...  

        # XXX: this is not necessarily correct. For example, it is *possible*
        # that a cluster might contain more than half of the genes and only
        # a few conditions.  Think about a "better" way to do this... i.e.
        # perhaps a matrix multiplication

        # returns which indices correspond to a cluster
        # TODO: fix this to get the cluster with the most non-zero entries
        cutKcluster <- function(mat)
        {
            cuts <- kmeans(mat, 2)$cluster

            cutIdx <- 1
            if (mean(mat[which(cuts == 1), ]) < 
                mean(mat[which(cuts == 2), ]))
            {
                cutIdx <- 2
            }
            cutIdx <- which(cuts == cutIdx)

            return(cutIdx)
        }

        rowIdx <- cutKcluster(A)
        colIdx <- cutKcluster(D)

        clusters <- append(clusters, list(list(rowIdx = rowIdx, colIdx = colIdx)))

        # given a list of rows and columns, converts pairwise combinations into
        # 1D index for matrix
        get1DIdx <- function(rows, cols) 
        {
            allPairs <- expand.grid(rows, cols)
            nrow(geneDf) * (allPairs[, 2] - 1) + allPairs[, 1]
        }

        # mask each gene with random values from the rest of the matrix also
        # considered doing this with JUST the other values of the gene... need
        # to experiment with that
        clustIdx <- get1DIdx(rowIdx, colIdx)
        # geneDf[clustIdx] <- sample(geneDf[-clustIdx], length(clustIdx))
        # temporarily sample rnorm to debug FP issue
        geneDf[clustIdx] <- rnorm(length(clustIdx))
    }

    return(list(allBcSol = allBcSol, clusters = clusters))
}

plotClusterExpression <- function(dat, clust)
{
    ggPlotExpression(dat[clust$rowIdx, clust$colIdx])
}



#' XXX: OBSOLETE! This is mostly still here for testing... 
biclustering <- function(geneDf, nSamples = 100)
{
    it <- 1
    aVal <- list()
    bVal <- list()
    dVal <- list()
    abVal <- list()
    splitIdx <- list()
    while (it <= nSamples) 
    {
        cat("Biclustering iteration: ", it, "\n")
        curSol <- maximizeOneSplit(geneDf)
        aVal <- append(aVal, list(curSol$a))
        bVal <- append(bVal, list(curSol$b))
        dVal <- append(dVal, list(curSol$d))
        abVal <- append(abVal, list(curSol$ab))
        splitIdx <- append(splitIdx, list(curSol$splitIdx))
        it <- it + 1
    }

    return(list(a = aVal, b = bVal, ab = abVal, d = dVal, splitIdx = splitIdx))
}

# load gene data supplied by Ben
# data is loaded in gene.mat
source("geneData.R", echo = FALSE)

# data from geneData.R
hi <- splitEvenly(nrow(gas.mat))
X.mat <- gas.mat[hi[[1]],] + 1
Y.mat <- gas.mat[hi[[2]],] + 1
X.mat <- log(X.mat)
Y.mat <- log(Y.mat)

X.mat.2 <- X.mat[1:20, 1:10]
Y.mat.2 <- Y.mat[1:20, 1:10]

startVec <- runif(ncol(X.mat.2))

meow2.optim <- biClustMax.optim(X.mat.2, Y.mat.2, 
                                startVec,
                                runif(nrow(X.mat.2), min = -1, max = 1), 
                                runif(nrow(Y.mat.2), min = -1, max = 1),
                                6
                                )
meow2.optim$d$par


# testing cgalC
Dmat <- matrix(c(2, 0, 0, 8), nrow = 2)
Amat <- matrix(c(1, 1, -1, 2), nrow = 2, byrow = T)
bvec <- c(7, 4)
lwr <- c(0, 0)
upr <- c(10, 4)
cvec <- c(0, -32)
c0 <- 64
cgalC(Dmat, Amat, bvec, lwr, upr, cvec, c0)



# attempting to implement using the nlminb function 
# d is a vector
# q is a vector
objectiveQ <- function(d, qVec, lam)
{
    sum(0.5 * d^2 * qVec + lam * sum(abs(d)))
}

startVec <- runif(nrow(meow2))

lapply(seq(0.05, 3, by = 0.05), function(lam) {
       nlminb(startVec, objectiveQ, 
              control = list(step.min = 0.05, step.max = 0.10),
              lower = rep.int(0, length(startVec)), 
              upper = rep.int(1, length(startVec)),
              qVec = diag(meow2), lam = -lam)$par
                           })


nlminb(startVec, objectiveQ, 
       control = list(step.min = 0.05, step.max = 0.10),
       lower = rep.int(0, length(startVec)), 
       upper = rep.int(1, length(startVec)),
       qVec = diag(meow2), lam = 1)

meow2vec <- diag(meow2)
pos <- rep(0, length(meow2vec))
pos[which(sign(meow2vec) == -1)] <- 1

objectiveQ(pos, meow2vec, -10)

# testing constrOptim
testConstrOptim <- function(q, lam, startPos)
{
    constrMat<- diag(length(q))
    constrMat <- rbind(constrMat, -diag(length(q)))
    constrMat <- rbind(constrMat, rep.int(-1, length(q)))
    constrLimit <- c(rep(0, length(q)), 
                     rep(-1, length(q)),
                     -lam)
    objectiveFn <- function(d, qVec)
    {
        sum(0.5 * d^2 * qVec)
    }
    constrOptim(startPos, objectiveFn, grad = NULL,
                ui = constrMat, ci = constrLimit,
                qVec = q
                )
}

testConstrOptim(meow2vec, 8, startVec)$par


# driver
sol <- maximizeOneSplit(rbind(X.mat.2, Y.mat.2))
bcSol <- biclustering(rbind(X.mat.2, Y.mat.2), nSamples = 10)

# FIXME: Fix the way it returns... basically make sure return is correct
