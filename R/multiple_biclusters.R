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
