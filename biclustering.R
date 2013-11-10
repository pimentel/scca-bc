# Author: Harold Pimentel
# Contact: pimentel@cs.berkeley.edu

# Lee et. al. package. Provides 'scca()'
library(scca)
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
                    penalty = "LASSO",
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
    constrLimit <- c(rep(0, length(q)), 
                     rep(-1, length(q)),
                     -lam)
    objectiveFn <- function(d, qVec)
    {
        # sum(0.5 * d^2 * qVec)
        sum(d^2 * qVec)
    }
    optimRes <- constrOptim(d.start, objectiveFn, grad = NULL,
                            control = list(maxit = 10000),
                            ui = constrMat, ci = constrLimit,
                            qVec = q
                            )

    if (optimRes$convergence != 0)
    {
        warning("Error with convergence.\n", "\tError code: ", 
                optimRes$convergence, "\tMessage: ", optimRes$message,
                immediate = TRUE)
    }

    return(list(a = a, b = b, d = optimRes$par, q))
}


##' Driver function to find solution for 
maximizeOneSplit <- function(geneDf, lam, epsA = 0.1, epsB = 0.1, epsD = 0.1, maxIt = 100)
{
    splitIdx <- splitEvenly(nrow(geneDf))
    xDf <- geneDf[splitIdx[[1]], ]
    yDf <- geneDf[splitIdx[[2]], ]

    # TODO: implement cross validation for lambda
    # lam <- round(min(dim(xDf)) * .3)
    # lam <- 20 * 2

    # randomly initialize d
    d <- 0
    maxUnif <- 1
    repeat {
        d <- runif(ncol(xDf), max = maxUnif)
        if (sum(d) <= lam)
            break
        else
            maxUnif <- min(0.05, maxUnif - 0.05)
    }
    # NB: values of a and b are not important currently since evaluated after d is set.
    # If ever change the order, fix this.
    a <- b <- runif(nrow(xDf), min = -1, max = 1)
    it <- 1
    curSol <- 0
    repeat {
        cat("\tOne split iteration: ", it, "\n")
        curSol <- biClustMax.optim(xDf, yDf, d, a, b, lam)
        cat ("\t\tdist: ",
             dist(rbind(curSol$a, a))[1], "\t",
             dist(rbind(curSol$b, b))[1], "\t",
             dist(rbind(curSol$d, d))[1], "\n")
        if (dist(rbind(curSol$a, a))[1] < epsA &
            dist(rbind(curSol$b, b))[1] < epsB &
            dist(rbind(curSol$d, d))[1] < epsD)
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
                splitIdx = splitIdx))
}

biclusteringPar <- function(geneDf, nSamples = 100, lam)
{
    mclapply(1:nSamples, function(it) {
             cat("Biclustering iteration: ", it, "\n")
             curSol <- maximizeOneSplit(geneDf, lam)
             return(curSol)
            #}, mc.silent = FALSE)
            })
}


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
