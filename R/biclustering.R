#' Randomly split a set of indices
#'
#' Randomly splits indices if there are an even number. If there is an odd
#' number, it will stop().
#' 
#' @param n the number of indices
#' @return a list with two disjoint sets of indices from 1 to n.
split_evenly <- function(n)
{
    if ( n %% 2 != 0)
    {
        stop("Not an even number of samples!")
    }
    sets <- split(sample.int(n), 1:2)
    return(sets)
}

#' Optimizing the biclustering function
#'
#'
#' One iteration of biclustering maximization. Takes the current values of 
#' d, and if successful returns one iteration.
#' 
#' @param X matrix with rows representing genes and columns representing conditions (n x k)
#' @param Y matrix with rows representing genes and columns representing conditions (n x k)
#' @param d.start vector of dimension k
#' @param lam maximum cluster size
#' @return A list containing maximum values of a, b, d
#' @export
biClustMax.optim <- function(X, Y, d.start, lam, verbose = TRUE, lam.lwr = 3.5,
    clustOptions = list())
{
    lamx <- 1:3
    if (!is.null(clustOptions$lamx))
        lamx <- clustOptions$lamx

    # First maximize a and b using SCCA
    res <- features_max_fscca(X, Y, d.start, lamx, lamx)
    a <- res$a
    b <- res$b

    d_optim <- lasso_max_d(X, Y, res$a, res$b, lam)


    return(list(a = a, b = b, d = d_optim, q, sccaLam = res$lambda))
}

##' Driver function to find solution for 
maximizeOneSplit <- function(exp_mat, lam, epsA = 0.001, epsB = 0.001, epsD = 0.01, maxIt = 100, lam.lwr, 
                             clustOptions = list())
{
    splitIdx <- split_evenly(nrow(exp_mat))
    X <- exp_mat[splitIdx[[1]], ]
    Y <- exp_mat[splitIdx[[2]], ]

    # randomly initialize d
    # FIXME: If there are lots of conditions, possible that won't be able to
    # find d that satisfy the constraint (if lam is sufficiently small. Come up
    # with different way to randomly assign values
    d <- 0
    minUnif <- 0.05
    maxUnif <- 1

    repeat { 
        d <- runif(ncol(X), max = maxUnif)
        if (sum(d) <= lam)
            break
        else
            maxUnif <- min(minUnif + 0.01, maxUnif - 0.05)
    }

    # NB: values of a and b are not important currently since evaluated after d is set.
    # If ever change the order, fix this.
    a <- b <- runif(nrow(X), min = -1, max = 1)
    it <- 1
    curSol <- list()
    cat(sprintf("% 10s% 7s% 7s\n", "a", "b", "d"))
    repeat {
        # cat("\tOne split iteration: ", it, "\n")
        curSol <- biClustMax.optim(X, Y, d, lam, lam.lwr = lam.lwr, clustOptions = clustOptions)

        a_dist <- dist(rbind(curSol$a, a))[1]
        b_dist <- dist(rbind(curSol$b, b))[1]
        # d_dist <- dist(rbind(curSol$d, d))[1]
        d_dist <- mean_absolute_tol(d, curSol$d)

        cat(sprintf("\t%.4f\t%.4f\t%.4f\n", a_dist, b_dist, d_dist))

        if (a_dist < epsA && b_dist < epsB && d_dist < epsD && it >= 3) 
        {
            cat("\tConverged. Number of iterations: ", it, "\n")
            break
        }

        a <- curSol$a
        b <- curSol$b
        d <- curSol$d
        it <- it + 1
        if (it > maxIt)
        {
            cat("Warning: exceeded maximum number of iterations (", maxIt, ")\n")
            break
        }
    }

    # NB: If there is something funky with the distribution of a and b, 
    # check here... though this should work correctly
    ab <- rep.int(0, length(b)*2)
    ab[splitIdx[[1]]] <- curSol$a
    ab[splitIdx[[2]]] <- curSol$b

    return(list(a = a, b = b, ab = ab, d = round(curSol$d, digits = 4),
                splitIdx = splitIdx, sccaLam = curSol$sccaLam))
}


#' Serial biclustering
#'
#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(cores = N).
#'
#' @param exp_mat matrix with features on rows and conditions on columns
#' @param nSamples integer denoting the number of permutations to perform
#' @param lam regularization parameter for conditions (maximum number of conditions allows in a bicluster)
#' @export
biclusteringPar <- function(exp_mat, nSamples = 100, lam, lam.lwr = 3.5, clustOptions = list())
{
    if (!is.matrix(exp_mat))
        stop("biclustering requires a matrix")

    mclapply(1:nSamples, function(it) {
             cat("Biclustering iteration: ", it, "\n")
             curSol <- maximizeOneSplit(exp_mat, lam = lam, lam.lwr = lam.lwr, clustOptions = clustOptions)
             return(curSol)
            })
}

#' Serial biclustering
#'
#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(cores = N).
#'
#' @param exp_mat matrix with features on rows and conditions on columns
#' @param nSamples integer denoting the number of permutations to perform
#' @param lam regularization parameter for conditions (maximum number of conditions allows in a bicluster)
#' @export
biclusteringSerial <- function(exp_mat, nSamples = 100, lam, lam.lwr = 3.5, clustOptions = list())
{
    if (!is.matrix(exp_mat))
        stop("biclustering requires a matrix")

    lapply(1:nSamples, function(it) {
        cat("Biclustering iteration: ", it, "\n")
        curSol <- maximizeOneSplit(exp_mat, lam = lam, lam.lwr = lam.lwr, clustOptions = clustOptions)
        return(curSol)
            })
}


#' Serial biclustering
#'
#' Uses one core to perform SCCAB
#'
#' @param exp_mat matrix with features on rows and conditions on columns
#' @param nSamples integer denoting the number of permutations to perform
#' @param lam regularization parameter for conditions (maximum number of conditions allows in a bicluster)
#' @export
bcSubSampleSerial <- function(exp_mat, nSamples = 100, lam, propSample = 0.6, lam.lwr = 3.5,
                           clustOptions = list())
{
    if (!is.matrix(exp_mat))
        stop("biclustering requires a matrix")

    if (propSample > 1 | propSample < 0.01)
        stop("Invalid range for propSample")

    nRowsSample <- round(propSample * nrow(exp_mat) )
    if (nRowsSample %% 2)
        nRowsSample <- nRowsSample + 1

    cat("Sampling ", nRowsSample, " features\n")
    lapply(1:nSamples, function(it) {
           cat("**Biclustering iteration: ", it, "\n")
           sampIdx <- sample.int(nrow(exp_mat), size = nRowsSample)
           abSol <- rep.int(NA, nrow(exp_mat))
           curSol <- maximizeOneSplit(exp_mat[sampIdx,], lam, lam.lwr = lam.lwr, 
                                      clustOptions = clustOptions)
           abSol[sampIdx] <- curSol$ab
           d <- curSol$d
           return(list(ab = abSol, d = d, sccaLam = curSol$sccaLam))
                           })
}


#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(cores = N).
#'
#' @param exp_mat matrix with genes on rows and columns defining conditions
#' @param nSamples integer denoting the number of permutations to perform
#' @param lam regularization parameter for conditions (maximum number of conditions allows in a bicluster)
#' @export
bcSubSamplePar <- function(exp_mat, nSamples = 100, lam, propSample = 0.6, lam.lwr = 3.5,
                           clustOptions = list())
{
    if (propSample > 1 | propSample < 0.01)
        stop("Invalid range for propSample")

    nRowsSample <- round(propSample * nrow(exp_mat) )
    if (nRowsSample %% 2)
        nRowsSample <- nRowsSample + 1

    cat("Sampling ", nRowsSample, " features\n")
    mclapply(1:nSamples, function(it) {
             cat("**Biclustering iteration: ", it, "\n")
             sampIdx <- sample.int(nrow(exp_mat), size = nRowsSample)
             abSol <- rep.int(NA, nrow(exp_mat))
             curSol <- maximizeOneSplit(exp_mat[sampIdx,], lam, lam.lwr = lam.lwr, 
                                        clustOptions = clustOptions)
             abSol[sampIdx] <- curSol$ab
             d <- curSol$d
             return(list(ab = abSol, d = d, sccaLam = curSol$sccaLam))
                           })
}
