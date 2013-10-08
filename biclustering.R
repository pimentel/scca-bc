# Author: Harold Pimentel
# Contact: pimentel@cs.berkeley.edu

# Provides 'ipop()' for quadratic programming
library(kernlab)
# Lee et. al. package. Provides 'scca()'
library(scca)
# Provides 'quadprog()'
library(quadprog)

##' One iteration of biclustering maximization. Takes the current values of 
##' d, a, b and if successful returns one iteration
##' 
##' @param xDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param yDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param d vector of dimension k
##' @param a vector of dimension n
##' @param b vector of dimension n
##' @return A list containing maximum values of a, b, d
biClustMax <- function(xDf, yDf, d, a, b, maxIt = 5, verbose = TRUE)
{
    quadRes <- 0
    quadIt <- 0

    # iterate until we've hit maximum iterations or QP has converged
    while (quadIt < maxIt & class(quadRes) == "numeric")
    {
        quadIt <- quadIt + 1
        if (verbose)
            cat("Computing tilde matrices\n")
        xTilde <- t(xDf %*% diag(d))
        yTilde <- t(yDf %*% diag(d))

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
        #
        # FIXME: Q ends up being singular. Way to solve is by making d sparse
        Q <- - 2 * diag(as.numeric((a %*% xDf) * (b %*% yDf)))
        return(Q)
        if (verbose)
        {
            cat("Solving quadratic program\n")
        }

        exitStatus <- tryCatch(
                               quadRes <- ipop(c = rep(0, length(d)),
                                               H = Q,
                                               A = diag(length(d)),
                                               b = rep(0, length(d)),
                                               l = rep(0, length(d)),
                                               u = rep(1, length(d)),
                                               r = rep(1, length(d))
                                               ),
                               error = function(e) {
                                   message(paste("ipop message:\n\t", e))
                                   return(1)
                               },
                               finally = {
                                   message(paste("Quadratic program completed iteration", quadIt, "/", maxIt))
                               })
    }

    if (exitStatus)
        stop("Non-zero exit status from 'kernlab.' Results from SCCA make 'H' matrix singular (see 'ipop' call).")
    if (verbose)
        cat("Successfully found solutation to QP\n")

    return(list(a = a, b = b, d = quadRes@primal))
}


##' One iteration of biclustering maximization. Takes the current values of 
##' d, a, b and if successful returns one iteration
##' 
##' @param xDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param yDf matrix with rows representing genes and columns representing conditions (n x k)
##' @param d vector of dimension k
##' @param a vector of dimension n
##' @param b vector of dimension n
##' @return A list containing maximum values of a, b, d
biClustMax.nlminb <- function(xDf, yDf, d, a, b, maxIt = 5, verbose = TRUE)
{

    # iterate until we've hit maximum iterations or QP has converged
    if (verbose)
        cat("Computing tilde matrices\n")
    xTilde <- t(xDf %*% diag(d))
    yTilde <- t(yDf %*% diag(d))

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

    objectiveQ <- function(d, lam = 0)
    {
        sum(0.5 * d^2 * q - lam * sum(abs(d)))
    }

    startVec <- runif(length(q))
    optimRes <- nlminb(startVec, objectiveQ, 
                       lower = rep.int(0, length(startVec)), 
                       upper = rep.int(1, length(startVec)),
                       lam = 0)

    return(list(a = a, b = b, d = optimRes, q))
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
        warning("Error with convergence. Check\n", "\tError code: ", 
                optimRes$convergence, "\tMessage: ", optimRes$message)
    }

    return(list(a = a, b = b, d = optimRes$par, q))
}


quadraticMax <- function(QMat, c, bound)
{
    dDim <- dim(QMat)[1]
    Amat <- matrix(0, nrow = dDim, ncol = dDim)
    sol <- ipop(c = rep(c, dDim),
                H = QMat,
                A = Amat,
                b = rep(0, dDim),
                l = rep(0, dDim),
                u = rep(1000, dDim),
                r = rep(0, dDim),
                bound = bound
                #margin = 0.05
                )

    return(sol)
}


# same as quadraticMax, but puts constraints on cluster size.
quadraticMax.clusterSize <- function(QMat, c, bound, clusterSize)
{
    dDim <- dim(QMat)[1]
    constraintMat <- matrix(0, nrow = dDim, ncol = dDim)
    constraintMat[1,] <- 1
    sol <- ipop(c = rep(c, dDim),
                H = QMat,
                A = constraintMat,
                b = rep(0, dDim),
                r = c(clusterSize, rep(0, dDim - 1)),
                l = rep(0, dDim),
                u = rep(1, dDim),
                bound = bound
                )
    return(sol)
}


quadProgMax <- function(DMat, lam)
{
    dDim <- dim(DMat)[1]
    # constraints that values must be greater than zero
    AMat <- diag(dDim)
    b0 <- rep(0, dDim)
    # constraints that values must be less than 1
    AMat <- rbind(AMat, -diag(dDim))
    b0 <- c(b0, rep(-1, dDim))
    # constraint that the sum of them are less than lam (regularization)
    # AMat <- rbind(AMat, rep(1, dDim))
    # b0 <- c(b0, lam)
    # dimensions are weird...
    AMat <- t(AMat)

    sol <- solve.QP(DMat, rep(lam, dDim), AMat, b0)

    return(sol)
}

library(QP)

cgalMax <- function(DMat, lam)
{
    dDim <- dim(DMat)[1]
    AMat <- matrix(0, nrow = 1, ncol = dDim)
    # AMat[1,] <- 1
    #sol <- cgalC(DMat, AMat, c(lam, rep(0, dDim-1)), 
    sol <- cgalQP(DMat, AMat, c(0), 
          rep(0, dDim), rep(1, dDim),
          rep(lam, dDim), 0)

    return(sol)
}

debug(quadProgMax)
undebug(quadProgMax)

quadProgMax(diag(x = 4, 3), 5)

quadProgMax(meow2, 5)

cgalMax(meow2, 1)

cgalMax(diag(5), 1)


solve.QP(diag(3), 1:3, matrix(0, 3, 3), rep(0, 3))

hi <- quadraticMax(meow, 0, 10)

hi <- quadraticMax(diag(5), 0, 10)

hi <- cgalMax(diag(5), 3)

hi <- quadraticMax(meow, -100000, 0.5)

hi <- quadraticMax(meow2, 0, 10)
hi <- quadraticMax(meow, -100000, 0.5)

quadraticMax(meow, 10000000, 100)
quadraticMax(meow, abs(0.0001*sum(meow)), 100)

sol <- quadraticMax.clusterSize(meow, 0, 10, max(abs(meow)))

debug(biClustMax)
undebug(biClustMax)

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


biClustMax(X, Y, diag(Dx), 1, 1)

##' Driver function to find solution for 
maximizeOneSplit <- function(geneDf)
{
    ## TODO: write me! in particular, iterate and check all exceptions
    splitIdx <- splitEvenly(nrow(geneDf))
    xDf <- geneDf[splitIdx[[1]], ]
    yDf <- geneDf[splitIdx[[2]], ]
    d <- runif(ncol(xDf))

    biClustMax(xDf, yDf, d, 1, 1)
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
