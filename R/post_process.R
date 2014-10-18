# XXX: This is the best performing method
#' @export
postSubSample.pca <- function(subSampleSol, abThresh = 0.6, dQuant = 0.5)
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    bootStrapNAs <- function(dat, lwr, upr)
    {
        for (icol in 1:ncol(dat))
        {
            rng <- quantile(dat[,icol], probs = c(lwr, upr), na.rm = T)
            # FIXME: why is this not matching at least some things?
            whichValid <- which( (rng[1] <= dat[,icol]) & (dat[,icol] <= rng[2]))
            nNa <- sum(is.na(dat[,icol]))
            if (length(whichValid) == 0)
            {
                # cat("len == 0", icol, "\n")
                # print(sum(!is.na(dat[,icol])))
                whichValid <- which(!is.na(dat[,icol]))
                # print(whichValid)
                # print("----")
            }
            samps <- NA
            if (length(whichValid) == 1)
                samps <- rep.int(whichValid, nNa)
            else
                samps <- sample(whichValid, nNa, replace = T)
                
            dat[is.na(dat[,icol]), icol] <- dat[samps,icol]
            # if (any(is.na(dat[,icol])))
            # {
            #     cat("samples:", icol, "\n")
            #     print(dat[samps,icol])
            #     print("whichValid values:")
            #     print(dat[whichValid,icol])
            #     print("whichValid:")
            #     print(length(whichValid))
            #     print(whichValid)
            # }
        }

        return(dat)
    }

    if (sum(is.na(ab)) > 0)
        ab <- bootStrapNAs(ab, 0.05, 0.95)

    pcaAB <- prcomp(ab, center = F, scale. = F)
    # pcaAB <- prcomp(ab)
    pcaAB$rotation <- as.data.frame(data.matrix(as.data.frame(pcaAB$rotation, stringsAsFactors = F)))

    pcaAB$rotation[,"mean"] <- apply(ab, 2, mean, na.rm = T)
    pcaAB$rotation[,"sd"] <- apply(ab, 2, sd, na.rm = T)
    pcaAB$rotation[,"median"] <- apply(ab, 2, median, na.rm = T)
    pcaAB$rotation[,"var"] <- apply(ab, 2, var, na.rm = T)
    pcaAB$rotation[,"cov"] <- with(pcaAB$rotation, sd / mean)

    pcaD <- prcomp(d, center = F, scale. = F)
    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)


    # colIdx <- clustKmeans(cbind(dQuant, dSd))
    # rowIdx <- clustKmeans(with(pcaAB$rotation, cbind(-PC1 * sdev[1], PC2 * sdev[2])))

    df <- cbind(-pcaAB$rotation$PC1 * pcaAB$sdev[1], 
                                pcaAB$rotation$PC2 * pcaAB$sdev[2])
    rowIdx <- clustKmeans(df)
    df2 <- cbind(-pcaD$rotation[,"PC1"] * pcaD$sdev[1],
        rnorm(nrow(pcaD$rotation), sd = 0.0001))
                 # pcaD$rotation[,"PC2"] * pcaD$sdev[2])
    colIdx <- clustKmeans(df2)
    # colIdx <- clustKmeans(cbind(dQuant))

    pcaAB$rotation$cluster <- FALSE
    pcaAB$rotation$cluster[rowIdx] <- TRUE

    pcaD$rotation <- data.frame(pcaD$rotation)
    pcaD$rotation[,"cluster"] <- "background"
    pcaD$rotation[colIdx,"cluster"] <- "bicluster"

    return(list(abDat = data.frame(pcaAB$rotation, sdev = pcaAB$sdev), dDat = pcaD,
            rowIdx = rowIdx, colIdx = colIdx))
    # return(list(rowIdx = rowIdx, colIdx = colIdx))
}

clustKmeans <- function(dat, minK = 2)
{
    gap <- cluster::clusGap(as.matrix(dat), kmeans, 5)
    nk <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
    if (nk == 1)
        nk <- minK
    cat("Conditions k: ", nk, "\n")
    kRes <- kmeans(dat, nk, nstart = 20)
    kMax <- which.max(kRes$centers[,1])

    which(kRes$cluster == kMax)
}

#' @export
postSubSample.ranks <- function(subSampleSol)
{
    ab <- getA(subSampleSol)
    ab_rank <- apply(ab, 2, rank, na.last = "keep")
    print(dim(ab))

    pca_ab <- prcomp(t(ab_rank), center = F, scale. = F)
    ab_features <- cbind(-pca_ab$rotation[,1] * pca_ab$sdev[1],
        pca_ab$rotation[,2] * pca_ab$sdev[2])
    rowIdx <- clustKmeans(ab_features)
    ab_features <- data.frame(ab_features)
    colnames(ab_features) <- c("PC1", "PC2")
    ab_features$in_clust <- FALSE
    ab_features$in_clust[rowIdx] <- TRUE


    d <- getD(subSampleSol)
    d_rank <- apply(d, 2, rank)
    pca_d <- prcomp(t(d), center = F, scale. = F)

    d_features <- cbind(-pca_d$rotation[,1] * pca_d$sdev[1],
        pca_d$rotation[,2] * pca_d$sdev[2])
    colIdx <- clustKmeans(d_features)
    d_features <- data.frame(d_features)
    colnames(d_features) <- c("PC1", "PC2")
    d_features$in_clust <- FALSE
    d_features$in_clust[colIdx] <- TRUE

    return(list(ab_rank = ab_rank, rowIdx = rowIdx, ab_features = ab_features,
            colIdx = colIdx, d_features = d_features))
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
        gap <- cluster::clusGap(as.matrix(dat), kmeans, 5)
        nk <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- minK
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }

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

    # rowIdx <- with(abStats, tight(cbind(medRate, sd)))
    # rowIdx <- with(abStats, tight(cbind(medRate)))
    # rowIdx <- with(abStats, tight(cbind(med, sd)))
    # rowIdx <- with(abStats, tight(cbind(med, sd)))

    clustKmeans <- function(dat, minK = 2)
    {
        gap <- cluster::clusGap(as.matrix(dat), kmeans, 5)
        nk <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
        if (nk == 1)
            nk <- minK
        cat("Conditions k: ", nk, "\n")
        kRes <- kmeans(dat, nk, nstart = 20)
        kMax <- which.max(kRes$centers[,1])

        which(kRes$cluster == kMax)
    }

    # colIdx <- clustKmeans(cbind(dQuant, dSd))
    colIdx <- clustKmeans(cbind(dQuant))

    return(list(abDat = abStats, cluster = list(rowIdx = rowIdx, colIdx = colIdx)))
}



postSubSample.tightPCA <- function(subSampleSol, abThresh = 0.6, dQuant = 0.5)
{
    # rows are different permutations,
    # columns are the genes or conditions
    ab <- t(sapply(subSampleSol, function (x) abs(x$ab)))
    d <- t(sapply(subSampleSol, function (x) x$d))

    bootStrapNAs <- function(dat, lwr, upr)
    {
        for (col in 1:ncol(dat))
        {
            rng <- quantile(dat[,col], probs = c(lwr, upr), na.rm = T)
            whichValid <- which(rng[1] <= dat[,col] & dat[,col] <= rng[2])
            nNa <- sum(is.na(dat[,col]))
            dat[is.na(dat[,col]), col] <- dat[sample(whichValid, nNa, replace = T), col]
        }

        return(dat)
    }

    # replace NAs w/ a bootstrap
    ab <- bootStrapNAs(ab, 0.05, 0.95)

    pcaAB <- prcomp(ab, center = F, scale. = F)
    pcaAB$rotation <- as.data.frame(data.matrix(as.data.frame(pcaAB$rotation, stringsAsFactors = F)))
    pcaAB$rotation[,"median"] <- apply(ab, 2, median, na.rm = T)
    pcaAB$rotation[,"mean"] <- apply(ab, 2, mean, na.rm = T)
    pcaAB$rotation[,"var"] <- apply(ab, 2, var, na.rm = T)
    pcaAB$rotation[,"sd"] <- apply(ab, 2, sd, na.rm = T)
    pcaAB$rotation[,"cov"] <- pcaAB$rotation[,"sd"] / pcaAB$rotation[,"mean"]

    dQuant <- apply(d, 2, quantile, probs = dQuant, na.rm = T)
    dSd <- apply(d, 2, sd, na.rm = T)

    # pcaAB$rotation[,c("PC2")] <- 0
    pcaAB$rotation[,"zeros"] <- 0

    # df <- pcaAB$rotation[,c("PC1", "PC2", "cov")]
    df <- pcaAB$rotation[,c("PC1", "zeros")]
    # df <- pcaAB$rotation[,c("PC1", "cov")]
    # df <- pcaAB$rotation[,c("PC1", "PC2")]


    # what if you remove the things that are so far away in distn
    # toInc <- which(df[,1] <= quantile(df[,1], probs = 0.5) )
    toInc <- which(df[,1] <= quantile(df[,1], probs = 1) )
    tc <- tight.clust(df[toInc,], 5, 10, standardize.gene = F)
    # tc <- tight.clust(df[toInc,], 2, 7, standardize.gene = F)
    pcaAB$rotation[,"cluster"] <- -1
    pcaAB$rotation[toInc,"cluster"] <- tc$cluster
    pcaAB$rotation[,"cluster"] <- as.factor(pcaAB$rotation[,"cluster"])
    clustCenter <- pcaAB$rotation %.% group_by(cluster) %.% summarize(center = mean(PC1)) 
    minClust <- which.min(clustCenter$center)
    # minClust <- 1

    # rowIdx <- which(1  == pcaAB$rotation[,"cluster"])
    rowIdx <- which(minClust  == pcaAB$rotation[,"cluster"])


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

    return(list(abDat = pcaAB$rotation, tc = tc$cluster, cluster = list(rowIdx = rowIdx, colIdx = colIdx)))
    return(list(rowIdx = rowIdx, colIdx = colIdx))
}


# Each row is a gene, each column is a different iteration
getA <- function(bcSol, takeAbs = TRUE)
{
    if (takeAbs)
        return( abs(sapply(bcSol, function (x) x$ab)) )

    sapply(bcSol, function (x) x$ab)
}

# Each row is a condition, each column is a different iteration
getD <- function(bcSol)
{
    sapply(bcSol, function (x) x$d)
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
