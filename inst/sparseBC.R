matrixBC.BIC <- function (x, k, r, lambda, alpha = 0.2, beta = 0.2, nstart = 20, 
    Sigma.init = NULL, Delta.init = NULL) 
{
    x <- x - mean(x)

    solz <- lapply(lambda, function(lam)
                     {
                         mclustering <- matrixBC(x, k, r, lambda = lam, 
                                                 nstart = nstart, alpha = alpha, beta = beta, Sigma.init = Sigma.init, 
                                                 Delta.init = Delta.init)
                         curBIC <- CalculateBICMatrix(x, mclustering)
                         nz <- sum(mclustering$Mus != 0)
                         list(BIC = curBIC, nonzero = nz, clust = mclustering)
                     })
    BIC <- sapply(solz, function(x) x$BIC)
    nonzero <- lapply(solz, function(x) x$nonzero)
    clust <- lapply(solz, function(x) x$clust)

    return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
        nonzeromus = nonzero, clusters = clust))
}

bictime <- system.time(matSol <- matrixBC.BIC(sol$mat, 2, 2, seq(0, 40, by= 2)))

bictime <- system.time(matSol <- sparseBC::matrixBC.BIC(hi, 2, 2, seq(0, 40, by= 2)))
fixInNamespace("matrixBC.BIC", "sparseBC")



sup <- mclapply(list(hi, hi, hi, hi, hi, hi), function(m)
         {
             sparseBC::matrixBC.BIC(m, 2, 2, seq(0, 40, by = 2))
         })


sparseBC.sim1.matrix <- function(data)
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 bicSol <- sparseBC::matrixBC.BIC(sol$permutedMat, 2, 2, c(0, 10))
                 bestLam <- 1
                 if(bicSol$lambda == 10)
                 {
                     bestLam <- 2
                 }
                 curSol <- bicSol$clusters[[bestLam]]
                 possRows <- list(which(curSol$Cs == 1), 
                                  which(curSol$Cs == 2))
                 possCols <- list(which(curSol$Ds == 1), 
                                  which(curSol$Ds == 2))
                 possRows <- lapply(possRows, function(x) perToOrigOrder(x, sol$rowOrder))
                 possCols <- lapply(possCols, function(x) perToOrigOrder(x, sol$colOrder))
                 allComb <- expand.grid(row = 1:2, col = 1:2)
                 allAmrs <- mapply(function(x, y) {
                                   amrs.hp(constMean.truth, 
                                           list(list(rowIdx = possRows[[x]], colIdx = possCols[[y]])))
                                  }, allComb$row, allComb$col)
                 bestComb <- which.max(allAmrs)
                 rowIdx <- unlist(possRows[allComb$row[bestComb]])
                 colIdx <- unlist(possCols[allComb$col[bestComb]])
                 list(sol = bicSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}

system.time(sparseBC.sim1.matrix.sol <- sparseBC.sim1.matrix(constMean.sol))
save(sparseBC.sim1.matrix.sol, file = "sparseBC.sim1.matrix.sol.RData")

sparseBC.sim3 <- function(data)
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 bicSol <- sparseBC::matrixBC.BIC(sol$permutedMat, 2, 2, c(0, 10))
                 bestLam <- 1
                 if(bicSol$lambda == 10)
                 {
                     bestLam <- 2
                 }
                 curSol <- bicSol$clusters[[bestLam]]
                 possRows <- list(which(curSol$Cs == 1), 
                                  which(curSol$Cs == 2))
                 possCols <- list(which(curSol$Ds == 1), 
                                  which(curSol$Ds == 2))
                 possRows <- lapply(possRows, function(x) perToOrigOrder(x, sol$rowOrder))
                 possCols <- lapply(possCols, function(x) perToOrigOrder(x, sol$colOrder))
                 allComb <- expand.grid(row = 1:2, col = 1:2)
                 allAmrs <- mapply(function(x, y) {
                                   amrs.hp(cor50Mvn.truth, 
                                           list(list(rowIdx = possRows[[x]], colIdx = possCols[[y]])))
                                  }, allComb$row, allComb$col)
                 bestComb <- which.max(allAmrs)
                 rowIdx <- unlist(possRows[allComb$row[bestComb]])
                 colIdx <- unlist(possCols[allComb$col[bestComb]])
                 list(sol = bicSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}

meow <- list(data = list(list(permutedMat = hi), list(permutedMat = hi), list(permutedMat = hi)))
# sparseBC.sim3(meow)
sim3.mat.time <- system.time(sparseBC.sim3.matrix <- sparseBC.sim3(sim50Mvn.sol))
save(sparseBC.sim3.matrix, file = "sparseBC.sim3.matrix.RData")

sim3.mat.time.highCor <- system.time(sparseBC.sim3.matrix.highCor.sol <- sparseBC.sim3(sim50Mvn.sol.highCor))
save(sparseBC.sim3.matrix.highCor.sol, file = "sparseBC.sim3.matrix.highCor.RData")

load("sparseBC.sim3.matrix.RData")
boxplot(computePerf(cor50Mvn.truth, sparseBC.sim3.matrix, amrs.hp), ylim = c(0,1))


