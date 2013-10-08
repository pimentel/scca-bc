
##' Function to generate normal noise and Gaussian blocks
##' 
##' @param nrow integer denoting the number of rows in the matrix
##' @param ncol integer denoting the number of columns in the matrix
##' @param noiseMean a mean for the background data
##' @param noiseSD a standard deviation for the background data
##' @param clusterOptions a list of lists. See examples.
##' @return A matrix of size nrow x ncol with blocks in locations denoted by clusterOptions
generateNormal <- function(nrow = 200, ncol = 40, noiseMean = 0, noiseSD = 4, 
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
        mat[clust$x.start:clust$x.end, clust$y.start:clust$y.end] <- 
            rnorm(clustNCols * clustNRows, mean = 20, sd = 4)
    }

    return(mat)
}

debug(generateNormal)
undebug(generateNormal)

testMat <- generateNormal(clusterOptions = 
                          list(list(x.start = 3, x.end = 10, y.start = 5, y.end = 20)))

testMat <- generateNormal(nrow = 40, ncol = 14, 
                          clusterOptions = 
                          list(list(x.start = 3, x.end = 10, y.start = 2, y.end = 3),
                               list(x.start = 13, x.end = 17, y.start = 1, y.end = 1)))
levelplot(testMat, col.regions=redgreen(75))
