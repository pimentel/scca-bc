library(s4vd)


load("constMean.sol.RData")

perToOrigOrder <- function(vars, permuted)
{
    permuted[vars]
}

test <- constMean.sol$data[[1]]$permutedMat
test2 <- test
test2[constMean.sol$data[[1]]$rowOrder, constMean.sol$data[[1]]$colOrder] <- test
all.equal(test2, constMean.sol$data[[1]]$mat)

curMat <- constMean.sol$data[[1]]$mat

curMat <- constMean.sol$data[[1]]$permutedMat
curSol <- biclust(curMat, method = BCPlaid(), cluster = "b", fit.model = y ~ m,
                  background.layer = TRUE)
curSol.bc <- list(list(rowIdx =  which(curSol@RowxNumber[,1]), colIdx = which(curSol@NumberxCol[1,])))

sort(perToOrigOrder(curSol.bc[[1]]$rowIdx, constMean.sol$data[[1]]$rowOrder))

sort(perToOrigOrder(curSol.bc[[1]]$colIdx, constMean.sol$data[[1]]$colOrder))

plaidComputeConst <- function(data) 
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 curSol <- biclust(sol$permutedMat, method = BCPlaid(), 
                                   cluster = "b", fit.model = y ~ m, background.layer = TRUE,
                                   max.layers = 1)
                 rowIdx <- which(curSol@RowxNumber[,1])
                 colIdx <- which(curSol@NumberxCol[1,])
                 rowIdx <- sort(perToOrigOrder(rowIdx, sol$rowOrder))
                 colIdx <- sort(perToOrigOrder(colIdx, sol$colOrder))
                 list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}


plaidComputeTwo <- function(data) 
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 curSol <- biclust(sol$permutedMat, method = BCPlaid(), 
                                   cluster = "b", fit.model = y ~ m + a + b, background.layer = TRUE,
                                   max.layers = 2)
                 if (curSol@Number == 2)
                 {
                     rowIdx1 <- which(curSol@RowxNumber[,1])
                     colIdx1 <- which(curSol@NumberxCol[1,])
                     rowIdx1 <- sort(perToOrigOrder(rowIdx1, sol$rowOrder))
                     colIdx1 <- sort(perToOrigOrder(colIdx1, sol$colOrder))
                     rowIdx2 <- which(curSol@RowxNumber[,2])
                     colIdx2 <- which(curSol@NumberxCol[2,])
                     rowIdx2 <- sort(perToOrigOrder(rowIdx2, sol$rowOrder))
                     colIdx2 <- sort(perToOrigOrder(colIdx2, sol$colOrder))
                     return(list(sol = curSol, clusters = list(list(rowIdx = rowIdx1, colIdx = colIdx1),
                                                        list(rowIdx = rowIdx2, colIdx = colIdx2))))
                 }
                 else
                 {
                     rowIdx <- which(curSol@RowxNumber[,1])
                     colIdx <- which(curSol@NumberxCol[1,])
                     rowIdx <- sort(perToOrigOrder(rowIdx, sol$rowOrder))
                     colIdx <- sort(perToOrigOrder(colIdx, sol$colOrder))
                     return(list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx))))
                 }
             })
}



constMean.truth <- list(list(rowIdx = 1:30, colIdx = 1:15))

system.time(plaid.constMean.sol <- plaidComputeConst(constMean.sol))
computePerf(constMean.truth, plaid.constMean.sol, amrs.hp)

system.time(ssvd.constMean.sol <- ssvdComputeConst(constMean.sol))
computePerf(constMean.truth, ssvd.constMean.sol, amrs.hp)
save(ssvd.constMean.sol, plaid.constMean.sol, file = "constMeanComp.RData")


load("sim50Mvn.sol.RData")
cor50Mvn.truth <- list(list(rowIdx = 1:300, colIdx = 1:30))
system.time(plaid.sim50Mvn <- plaidComputeConst(sim50Mvn.sol))
computePerf(cor50Mvn.truth, plaid.sim50Mvn, amrs.hp)

load("sim50Mvn.sol.highCor.RData")
cor50Mvn.truth <- list(list(rowIdx = 1:300, colIdx = 1:30))
system.time(plaid.sim50Mvn.highCor <- plaidComputeConst(sim50Mvn.sol.highCor))
system.time(ssvd.sim50Mvn.highCor <- ssvdComputeConst(sim50Mvn.sol.highCor))
save(plaid.sim50Mvn.highCor, ssvd.sim50Mvn.highCor, file = "mvn50comp.highCor.RData")

computePerf(cor50Mvn.truth, plaid.sim50Mvn, amrs.hp)

ssvdComputeConst <- function(data) 
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 curSol <- biclust(sol$permutedMat, method = BCs4vd(), 
                                   nbiclust = 1,
                                   row.overlap = FALSE, col.overlap = FALSE, iter = 500,
                                   cols.nc = T, rows.nc = T)
                 rowIdx <- which(curSol@RowxNumber[,1])
                 colIdx <- which(curSol@NumberxCol[1,])
                 rowIdx <- sort(perToOrigOrder(rowIdx, sol$rowOrder))
                 colIdx <- sort(perToOrigOrder(colIdx, sol$colOrder))

                 list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}


ssvdComputeConst2 <- function(data) 
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 curSol <- biclust(sol$permutedMat, method = BCs4vd(), 
                                   nbiclust = 2,
                                   row.overlap = FALSE, col.overlap = FALSE, iter = 500,
                                   cols.nc = T, rows.nc = T)

                 if (curSol@Number == 2)
                 {
                     rowIdx1 <- which(curSol@RowxNumber[,1])
                     colIdx1 <- which(curSol@NumberxCol[1,])
                     rowIdx1 <- sort(perToOrigOrder(rowIdx1, sol$rowOrder))
                     colIdx1 <- sort(perToOrigOrder(colIdx1, sol$colOrder))
                     rowIdx2 <- which(curSol@RowxNumber[,2])
                     colIdx2 <- which(curSol@NumberxCol[2,])
                     rowIdx2 <- sort(perToOrigOrder(rowIdx2, sol$rowOrder))
                     colIdx2 <- sort(perToOrigOrder(colIdx2, sol$colOrder))
                     return(list(sol = curSol, clusters = list(list(rowIdx = rowIdx1, colIdx = colIdx1),
                                                        list(rowIdx = rowIdx2, colIdx = colIdx2))))
                 }
                 else
                 {
                     rowIdx <- which(curSol@RowxNumber[,1])
                     colIdx <- which(curSol@NumberxCol[1,])
                     rowIdx <- sort(perToOrigOrder(rowIdx, sol$rowOrder))
                     colIdx <- sort(perToOrigOrder(colIdx, sol$colOrder))
                     return(list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx))))
                 }
             })
}




system.time(ssvd.constMean.sol <- ssvdComputeConst(constMean.sol))
computePerf(constMean.truth, ssvd.constMean.sol, amrs.hp)
computePerf(constMean.truth, ssvd.constMean.sol, amrs)

system.time(ssvd.sim50Mvn <- ssvdComputeConst(sim50Mvn.sol))
computePerf(cor50Mvn.truth, ssvd.sim50Mvn, amrs.hp)

save(ssvd.sim50Mvn, plaid.sim50Mvn, file = "mvn50comp.RData")


# deal with two BC now

load("twoBC.sol.RData")
truth.twoBC <- list(list(rowIdx = 1:50, colIdx = 1:15), 
                    list(rowIdx = 200:239, colIdx = 26:40))


plaid.twoBC.sol <- plaidComputeTwo(twoBC.sol)
ssvd.twoBC.sol <- ssvdComputeConst2(twoBC.sol)
save(plaid.twoBC.sol, ssvd.twoBC.sol, file = "twoBC.sol.comp.RData")

computePerf(truth.twoBC, plaid.twoBC.sol, amrs.hp)
computePerf(truth.twoBC, twoBC.sol$sol, amrs.hp)
computePerf(truth.twoBC, ssvd.twoBC.sol, amrs.hp)

load("twoBC.sol.comp.RData")
computePerf(list(truth.twoBC[[1]]), getClusterX(ssvd.twoBC.sol, 1), amrs.hp)
computePerf(list(truth.twoBC[[2]]), getClusterX(ssvd.twoBC.sol, 2), amrs.hp)

getClusterX(twoBC.sol$sol, 1)

computePerf(list(truth.twoBC[[1]]), getClusterX(plaid.twoBC.sol, 1), amrs.hp)
computePerf(list(truth.twoBC[[2]]), getClusterX(plaid.twoBC.sol, 2), amrs.hp)

computePerf(truth.twoBC, ssvd.twoBC.sol, amrs.eachClust)
computePerf(truth.twoBC, plaid.twoBC.sol, amrs.eachClust)

# debugging
curDat <- twoBC.sol$data[[1]]$permutedMat
meowSol <- biclust(curDat, method = BCs4vd(), 
                                   nbiclust = 2,
                                   row.overlap = FALSE, col.overlap = FALSE, iter = 500,
                                   cols.nc = T, rows.nc = T)


sparseBC.sim1 <- function(data)
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 bicSol <- sparseBC.BIC(sol$permutedMat, 2, 2, seq(0, 40, by = 2))
                 curSol <- sparseBC(sol$permutedMat, 2, 2, bicSol$lambda)
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
                 list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}

sparseBC.sim2 <- function(solData)
{
    mclapply(1:length(solData), function(idx)
             {
                 cat("Currently on IDX: ", idx, "\n")
                 sol <- solData[[idx]]
                 bicSol <- sparseBC.BIC(sol$permutedMat, 3, 3, seq(0, 40, by = 2))
                 curSol <- sparseBC(sol$permutedMat, 3, 3, bicSol$lambda)
                 combinations <- combn(9, 2, simplify = F)
                 allCombs <- expand.grid(row = 1:3, col = 1:3)
                 bcIdx <- lapply(combinations, function(comb)
                                 {
                                     clust1 <- list(rowIdx = which(curSol$Cs == allCombs$row[comb[1]]),
                                                    colIdx = which(curSol$Ds == allCombs$col[comb[1]]))
                                     clust2 <- list(rowIdx = which(curSol$Cs == allCombs$row[comb[2]]),
                                                    colIdx = which(curSol$Ds == allCombs$col[comb[2]]))
                                     list(clust1, clust2)
                                 })
                 bcIdx <- lapply(bcIdx, function(bc)
                                 {
                                     clust1 <- list(rowIdx = perToOrigOrder(bc[[1]]$rowIdx, sol$rowOrder),
                                                    colIdx = perToOrigOrder(bc[[1]]$colIdx, sol$colOrder))
                                     clust2 <- list(rowIdx = perToOrigOrder(bc[[2]]$rowIdx, sol$rowOrder),
                                                    colIdx = perToOrigOrder(bc[[2]]$colIdx, sol$colOrder))
                                     list(clust1, clust2)
                                 })
                 allAmrs <- sapply(bcIdx, function(bc)
                                   {
                                       amrs.hp(truth.twoBC, bc)
                                   })
                 bestAmrs <- which.max(allAmrs)
                 list(sol = curSol, clusters = bcIdx[[bestAmrs]],
                      firstClust = allCombs[combinations[[bestAmrs]],],
                      secClust = allCombs[combinations[[bestAmrs]],]
                      )
             })
}


sparseBC.sim3.normal <- function(data)
{
    solData <- data$data
    mclapply(solData, function(sol)
             {
                 bicSol <- sparseBC.BIC(sol$permutedMat, 2, 2, seq(0, 40, by = 2))
                 curSol <- sparseBC(sol$permutedMat, 2, 2, bicSol$lambda)
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
                 list(sol = curSol, clusters = list(list(rowIdx = rowIdx, colIdx = colIdx)))
             })
}


sparseBC.sim2.matrix <- function(solData)
{
    mclapply(1:length(solData), function(idx)
             {
                 cat("Currently on IDX: ", idx, "\n")
                 sol <- solData[[idx]]
                 bicSol <- matrixBC.BIC(sol$permutedMat, 3, 3, c(0, 10))
                 curSol <- matrixBC(sol$permutedMat, 3, 3, bicSol$lambda, 0.2, 0.2)
                 combinations <- combn(9, 2, simplify = F)
                 allCombs <- expand.grid(row = 1:3, col = 1:3)
                 bcIdx <- lapply(combinations, function(comb)
                                 {
                                     clust1 <- list(rowIdx = which(curSol$Cs == allCombs$row[comb[1]]),
                                                    colIdx = which(curSol$Ds == allCombs$col[comb[1]]))
                                     clust2 <- list(rowIdx = which(curSol$Cs == allCombs$row[comb[2]]),
                                                    colIdx = which(curSol$Ds == allCombs$col[comb[2]]))
                                     list(clust1, clust2)
                                 })
                 bcIdx <- lapply(bcIdx, function(bc)
                                 {
                                     clust1 <- list(rowIdx = perToOrigOrder(bc[[1]]$rowIdx, sol$rowOrder),
                                                    colIdx = perToOrigOrder(bc[[1]]$colIdx, sol$colOrder))
                                     clust2 <- list(rowIdx = perToOrigOrder(bc[[2]]$rowIdx, sol$rowOrder),
                                                    colIdx = perToOrigOrder(bc[[2]]$colIdx, sol$colOrder))
                                     list(clust1, clust2)
                                 })
                 allAmrs <- sapply(bcIdx, function(bc)
                                   {
                                       amrs.hp(truth.twoBC, bc)
                                   })
                 bestAmrs <- which.max(allAmrs)
                 list(sol = curSol, clusters = bcIdx[[bestAmrs]],
                      firstClust = allCombs[combinations[[bestAmrs]],],
                      secClust = allCombs[combinations[[bestAmrs]],]
                      )
             })
}




debug(sparseBC.sim2)
sol <- twoBC.sol$data[[1]]
system.time(sparseBC.sim2 <- sparseBC.sim2(list(twoBC.sol$data[[1]])))
save(sparseBC.sim2, file = "sparseBC.sim2.RData")

system.time(sparseBC.sim2.matrix.sol <- sparseBC.sim2.matrix(twoBC.sol$data))
save(sparseBC.sim2.matrix.sol, file = "sparseBC.sim2.matrix.sol.RData")

getClusterX <- function(solList, idx)
{
    lapply(solList, function(curSol)
           {
               clust <- curSol$clusters
               maxIdx <- which.max(lapply(clust, function(x) amrs.hp(list(truth.twoBC[[idx]]), list(x))))
               return(list(clusters = list(clust[[maxIdx]])))
           })
}
computePerf(truth.twoBC, sparseBC.sim2.matrix.sol, amrs.hp)

sparseBC.sim2.clust1 <- getClusterX(sparseBC.sim2, 1)
computePerf(list(truth.twoBC[[1]]), sparseBC.sim2.clust1, amrs.hp)

sparseBC.sim2.clust2 <- getClusterX(sparseBC.sim2, 2)
computePerf(list(truth.twoBC[[2]]), sparseBC.sim2.clust2, amrs.hp)

system.time(sparseBC.sim2 <- sparseBC.sim2(twoBC.sol$data))
computePerf(truth.twoBC, sparseBC.sim2)

system.time(sparseBC.sim3.normal.sol <- sparseBC.sim3.normal(sim50Mvn.sol))
save(sparseBC.sim3.normal.sol, file = "sparseBC.sim3.normal.RData")


system.time(sparseBC.sim3b.normal.sol <- sparseBC.sim3.normal(sim50Mvn.sol.highCor))
save(sparseBC.sim3b.normal.sol, file = "sparseBC.sim3b.normal.RData")

computePerf(cor50Mvn.truth, sparseBC.sim3.normal.sol, amrs.hp)

# end debugging
sparseBC.BIC(curMat, 2, 2, seq(0, 40, length = 10))
sparseBC(curMat, 2, 2, 0)

# sim 1 sparseBC (witten)
debug(sparseBC.sim1)
sol <- constMean.sol$data[[1]]
system.time(sparseBC.constMean.sol <- sparseBC.sim1(constMean.sol))
computePerf(constMean.truth, sparseBC.constMean.sol, amrs.hp)
save(sparseBC.constMean.sol, file = "sparseBC.sim1.RData")

load("sparseBC.sim1.RData")




system.time(sparseBC.sim2.sol <- sparseBC.sim2(twoBC.sol))

# sim 3 witten
system.time(sparseBC.sim3.sol <- sparseBC.sim1(sim50Mvn.sol))
sol <- sim50Mvn.sol$data[[1]]
matSol <- matrixBC.BIC(sol$mat, 2, 2, seq(0, 40, by= 2))
matSol <- matrixBC.BIC(sol$mat, 2, 2, c(0, 10, 20, 30, 40))
matrixBC.BICl()

