# Author: Harold Pimentel
# Contact: pimentel@cs.berkeley.edu

# Lee et. al. package. Provides 'scca()'
library(scca)

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
biClustMax.optim <- function(xDf, yDf, d.start, a, b, lam, verbose = TRUE)
{

    # iterate until we've hit maximum iterations or QP has converged
    if (verbose)
        cat("Computing tilde matrices\n")
    xTilde <- t(xDf %*% diag(d.start))
    yTilde <- t(yDf %*% diag(d.start))

    if (verbose)
        cat("Computing SCCA component\n")
    # First maximize a and b using SCCA
    sccaRes <- scca(as.matrix(xTilde), as.matrix(yTilde), 
                    nc = 1, 
                    penalty = "HL",
                    center = TRUE, scale = FALSE)
    a <- as.numeric(sccaRes$A)
    b <- as.numeric(sccaRes$B)

    # Now, maximize D using quadratic programming
    q <- - 2 * as.numeric((a %*% xDf) * (b %*% yDf))

    if (verbose)
        cat("Solving optimization program\n")

    # lower bound is zero
    constrMat<- diag(length(q))
    # upper bound is one
    constrMat <- rbind(constrMat, -diag(length(q)))
    # regularize (sum(d) <= lam)
    constrMat <- rbind(constrMat, rep.int(-1, length(q)))
    constrLimit <- c(rep(0, length(q)), 
                     rep(-1, length(q)),
                     -lam)
    objectiveFn <- function(d, qVec)
    {
        sum(0.5 * d^2 * qVec)
    }
    optimRes <- constrOptim(d.start, objectiveFn, grad = NULL,
                            ui = constrMat, ci = constrLimit,
                            qVec = q
                            )

    if (optimRes$convergence != 0)
    {
        warning("Error with convergence.\n", "\tError code: ", 
                optimRes$convergence, "\tMessage: ", optimRes$message)
    }

    return(list(a = a, b = b, d = optimRes$par, q))
}


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


##' Driver function to find solution for 
maximizeOneSplit <- function(geneDf)
{
    ## TODO: write me! in particular, iterate and check all exceptions
    splitIdx <- splitEvenly(nrow(geneDf))
    xDf <- geneDf[splitIdx[[1]], ]
    yDf <- geneDf[splitIdx[[2]], ]

    lam <- round(min(xDf) * .3)

    d <- 0
    repeat {
        d <- runif(ncol(xDf))
        if (sum(d) <= lam)
            break
    }

    # NB: values of a and b are not important currently since evaluated after
    aVal <- list()
    bVal <- list()
    dVal <- list()
    repeat {
        biClustMax(xDf, yDf, d, 1, 1, lam)

    }
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

# testing one iteration of biclustering
meow <- biClustMax(X.mat, Y.mat, 
                   rep(1, ncol(gas.mat)),
                   runif(length(hi[[1]]), min = -1, max = 1), 
                   runif(length(hi[[1]]), min = -1, max = 1)
                   )

# testing one iteration of biclustering
meow2 <- biClustMax(X.mat.2, Y.mat.2, 
                   rep(1, ncol(X.mat.2)),
                   runif(nrow(X.mat.2), min = -1, max = 1), 
                   runif(nrow(Y.mat.2), min = -1, max = 1)
                   )

meow.nlminb <- biClustMax.nlminb(X.mat, Y.mat, 
                                 runif(ncol(X.mat)),
                                 runif(nrow(X.mat), min = -1, max = 1), 
                                 runif(nrow(Y.mat), min = -1, max = 1)
                                 )

meow2.nlminb <- biClustMax(X.mat.2, Y.mat.2, 
                           runif(ncol(X.mat.2)),
                           runif(nrow(X.mat.2), min = -1, max = 1), 
                           runif(nrow(Y.mat.2), min = -1, max = 1)
                           )

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
