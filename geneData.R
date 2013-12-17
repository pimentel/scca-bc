geneData <- read.csv("../data/Supplementary_Data_File_10_gene_expression_BPKM.csv",
                     header = T)

geneDataAdult <- geneData[, -grep("L3", colnames(geneData))]
geneDataAdult <- geneDataAdult[, order(colnames(geneDataAdult))]
gas <- geneDataAdult[, order(colnames(geneDataAdult))]
gas.mat <- as.matrix(gas)

removeAllZeroRows <- function(mat)
{
    whichZero <- apply(mat, 1, function(row) all(row < .Machine$double.eps))
    whichNonZero <-  ! whichZero
    return(whichNonZero)
}

gas.mat.nz <- gas.mat[removeAllZeroRows(gas.mat),]


trimPercent <- function(mat, valid, per)
{
    tmpIdx <- apply(mat, 1, function(row)
                    {
                        mean(row > 0) >= per
                    })
    mat2 <- mat[tmpIdx,]
    validAtFDR <- valid[(valid %in% rownames(mat2))]
    return(validAtFDR)
}

idxAtFDR <- function(fdrMat, validNames, fdr)
{
    which.fdr <- apply(fdrMat, 1, function(row) any(row <= fdr))
    which.fdr <- names(which.fdr)[which(which.fdr == TRUE)]
    validIdx <- which.fdr[which((which.fdr %in% validNames))]
    return(validIdx)
}

meow <- trimPercent(gas.mat, validIdx, 0.9)
meow <- trimPercent(gas.mat, rownames(gas.mat), 0.9)

meow <- trimPercent(gas.mat, rownames(gas.mat), 0.95)

meow <- trimPercent(gas.mat, rownames(gas.mat), 1)

length(trimPercent(gas.mat, idxAtFDR(geneDE.fdr, rownames(gas.mat), 0.3), 0.9))

fdr30nz90.idx <- trimPercent(gas.mat, idxAtFDR(geneDE.fdr, rownames(gas.mat), 0.3), 0.9)
fdr30nz90.dat <- gasMatZ[fdr30nz90.idx,]
fdr30nz90.bc <- bcSubSamplePar(gasMatZ[fdr30nz90.idx,],100, 20, 0.6)
save(fdr30nz90.bc, file = "fdr30nz90.bc.RData")

load("fdr30nz90.bc.RData", verbose = T)
tmpDat <- gas.mat[fdr30nz90.idx,]
fdr30nz90.pss <- postSubSample.percent(fdr30nz90.bc, 0.95, 0.5)
fdr30nz90.pss
levelplot(fdr30nz90.dat[fdr30nz90.pss$rowIdx, fdr30nz90.pss$d],
         scales=list(x=list(rot=90)), col.regions = redgreen(75))
write(rownames(fdr30nz90.dat)[fdr30nz90.pss$rowIdx], "../results/fdr30nz90_FG.txt")
write(rownames(fdr30nz90.dat), "../results/fdr30nz90_BG.txt")


fdr30nz90.pss <- postSubSample.percent(fdr30nz90.bc, 0.90, 0.5)
fdr30nz90.pss
levelplot(fdr30nz90.dat[fdr30nz90.pss$rowIdx, fdr30nz90.pss$d],
         scales=list(x=list(rot=90)), col.regions = redgreen(75))
write(rownames(fdr30nz90.dat)[fdr30nz90.pss$rowIdx], "../results/fdr30nz90_90p_FG.txt")
write(rownames(fdr30nz90.dat), "../results/fdr30nz90_BG.txt")

tmpDat <- gas.mat[rownames(fdr30nz90.dat),]

levelplot(tmpDat[fdr30nz90.pss$rowIdx, fdr30nz90.pss$d],
         scales=list(x=list(rot=90)), col.regions = redgreen(75))

tmpIdx <- apply(gas.mat.nz, 1, function(row)
                {
                    mean(row > 0) >= 0.9
                })

gas.mat.75 <- gas.mat.nz[tmpIdx,]
#nrow(gas.mat.50[validIdx,])
remove75 <- validIdx[(validIdx %in% rownames(gas.mat.75))]

remove75 <- validIdx[(validIdx %in% rownames(gas.mat.75))]

length(remove75)
gasMat.75.bc <- bcSubSamplePar(gas.mat.75[remove75, ], 100, 20, 0.6)
save(gasMat.75.bc, file = "gasMat.75.RData")

load("gasMat.75.RData")

pss.75 <- postSubSample.percent(gasMat.75.bc, 0.90, 0.5)
rm75 <- gas.mat.75[remove75,]
levelplot(rm75[pss.75$rowIdx, pss.75$d], scales=list(x=list(rot=90)), 
          col.regions = redgreen(75))

length(validIdx[(validIdx %in% rownames(gas.mat.50))])

gasMatZ <- t(scale(t(gas.mat.nz)))
system.time(gasMatZ.bc <- bcSubSamplePar(gasMatZ[validIdx,], 100, 20, 0.6))
save(gasMatZ.bc, gasMatZ, file = "gasMatZ.RData")

load("gasMatZ.RData")
postGasZ <- postSubSample.percent(gasMatZ.bc, 0.90, 0.5)
tmpDat <- gasMatZ[validIdx,]
levelplot(tmpDat[postGasZ$rowIdx, postGasZ$colIdx], 
          scales=list(x=list(rot=90)), col.regions = redgreen(75))

tmpDat <- gas.mat.nz[validIdx,]
levelplot(tmpDat[postGasZ$rowIdx, postGasZ$colIdx], 
          scales=list(x=list(rot=90)), col.regions = redgreen(75))
quantile(tmpDat[postGasZ$rowIdx, postGasZ$colIdx])

# d <- diag(runif(ncol(gas)))
# colnames(d) <- colnames(gas.mat)


logGA <- log(gas.mat + 1)
logGA <- logGA[validIdx,]
whichZero <- apply(logGA, 1, function(row) all(row == 0.0))
whichNonZero <- ! whichZero
logGA <- logGA[whichNonZero, ]
logGA <- logGA[validIdx[-1539],]

whichZero <- apply(gas.mat, 1, function(row) all(row == 0.0))
whichNonZero <- !whichZero
gas.trim <- gas.mat[whichNonZero,]
gas.trim <- gas.trim[validIdx[-1539],]

gasBC.raw <- bcMultipleClusters(gas.trim, 20, 1, 100)

gasBC.raw <- bcSubSamplePar(gas.trim, 100, 20, 0.6)

gasBC.log <- bcSubSamplePar(log(gas.trim + 1), 100, 20, 0.6)
save(gasBC.log, file = "gasBC.log.RData")

load("gasBC.log.RData")
postGasLog <- postSubSample.percent(gasBC.log, 0.75, 0.75)
postGasLog
tmpData <- t(scale(t(gas.trim[validIdx[-1539],])))
tmpData2 <- (gas.trim[validIdx[-1539],])
levelplot(tmpData[postGasLog$rowIdx, postGasLog$colIdx],
          scales=list(x=list(rot=90)))

levelplot(tmpData2[postGasLog$rowIdx, postGasLog$colIdx],
          scales=list(x=list(rot=90)))


write(rownames(gas.trim[postGasLog$rowIdx, postGasLog$colIdx]), "../results/fdr20_log_FG.txt")
write(rownames(gas.trim[validIdx[-1539],]), "../results/fdr20_log_BG.txt")

corMat <- cor(t(gas.trim[postGasLog$rowIdx, postGasLog$colIdx]))
quantile(gas.trim[postGasLog$rowIdx, postGasLog$colIdx])

gasTrimZ <- t(scale(t(gas.trim)))
levelplot(gasTrimZ[postGasLog$rowIdx, postGasLog$colIdx],
          scales=list(x=list(rot=90)))

hist(abs(corMat[upper.tri(corMat)]))


postGasBC.raw <- postSubSample(gasBC.raw)
postGasBC.rawClust <- post.hclust(postGasBC.raw, 4, 2)
postGasBC.rawClust



save(gasBC.raw, file = "gasBC.raw.RData")

load("gasBC.raw.RData")
postGasBC <- postSubSample.percent(gasBC.raw, 0.90, 0.5)
gas.trim[postGasBC$rowIdx, postGasBC$d]

postGasBC <- postSubSample.percent(gasBC.raw, 0.75, 0.75)
postGasBC

levelplot(sqrt(gas.trim[postGasBC$rowIdx, postGasBC$colIdx]),
          scales=list(x=list(rot=90)))

hist(sqrt(gas.trim[postGasBC$rowIdx, postGasBC$colIdx]), breaks = 50)
table(gas.trim[postGasBC$rowIdx, postGasBC$colIdx])

plotParSolution(gasBC.raw$allBcSol[[1]])
gasBC.raw$clusters[[1]]


gasBC.rawHL <- bcMultipleClusters(gas.trim, 20, 1, 100)

save(gasBC.rawHL, file = "gasBC.rawHL.RData")

load("gasBC.rawHL.RData")

plotParSolution(gasBC.rawHL$allBcSol[[1]])
gasBC.rawHL$clusters[[1]]


gasBC.rawSCAD <- bcMultipleClusters(gas.trim, 20, 1, 100)
save(gasBC.rawSCAD, file = "gasBC.rawSCAD.RData")

load("gasBC.rawSCAD.RData")

plotParSolution(gasBC.rawSCAD$allBcSol[[1]])
gasBC.rawSCAD$clusters[[1]]

gasBC.rawSOFT <- bcMultipleClusters(gas.trim, 20, 1, 100)
save(gasBC.rawSOFT, file = "gasBC.rawSOFT.RData")

gasBC <- bcMultipleClusters(logGA[validIdx[-1539],], 20, 1, 100)
save(gasBC, file = "gasBC.RData")
load("gasBC.RData")

plotParSolution(gasBC$allBcSol[[1]])
absA <- abs(sapply(gasBC$allBcSol[[1]], function(x) x$ab))

gasBC.corCondRem <- bcMultipleClusters(logGA[,-gasBC$clusters[[1]]$colIdx], 20, 1, 100)
save(gasBC.corCondRem, file = "gasBC.corCondRem.RData")

load("gasBC.corCondRem.RData")

plotParSolution(gasBC.corCondRem$allBcSol[[1]])
gasBC.corCondRem$clusters[[1]]


levelplot(absA)
heatmap(absA, Colv = NA)

gaSplits <- splitEvenly(nrow(logGA))
gaX <- logGA[gaSplits[[1]], ]
gaY <- logGA[gaSplits[[2]], ]

scale.gaX <- scale(gaX)


# getting genes for which there is a low FDR
geneDE <-
    read.csv("../data/Supplementary_table_9_treatments_differential_expression (1).csv", header = T)
geneDE.fdr <- geneDE[,grep("FDR", colnames(geneDE))]
rownames(geneDE.fdr) <- geneDE[,1]

which.20.fdr <- apply(geneDE.fdr, 1, function(row) any(row <= 0.20))
which.20.fdr <- names(which.20.fdr)[which(which.20.fdr == TRUE)]
which.20.fdr[which(!(which.20.fdr %in% rownames(gas.mat)))]
validIdx <- which.20.fdr[which((which.20.fdr %in% rownames(gas.mat)))]


# starting from scratch !
# gas.mat is unadulterated
gasMatNz <- removeAllZeroRows(gas.mat)
gasMatNz <- gas.mat[gasMatNz,]

# getting next BC
get1DIdx.ext <- function(geneDf, rows, cols) 
{
    allPairs <- expand.grid(rows, cols)
    nrow(geneDf) * (allPairs[, 2] - 1) + allPairs[, 1]
}

# FDR 30, 90% of conditions contain values > 0
fdr30p90.idx <- trimPercent(gasMatNz, idxAtFDR(geneDE.fdr, rownames(gasMatNz), 0.30), 0.9)
fdr30p90 <- gasMatNz[fdr30p90.idx,]

bc30.fdr30p90.sqrt <- bcSubSamplePar(sqrt(fdr30p90), 100, 30, 0.6)
bc25.fdr30p90.sqrt <- bcSubSamplePar(sqrt(fdr30p90), 100, 25, 0.6)
bc20.fdr30p90.sqrt <- bcSubSamplePar(sqrt(fdr30p90), 100, 20, 0.6)
save(
     bc30.fdr30p90.sqrt,
     bc25.fdr30p90.sqrt,
     bc20.fdr30p90.sqrt, file = "bcFdr30p90.RData")


fdr20.idx <- idxAtFDR(geneDE.fdr, rownames(gasMatNz), 0.20)
fdr20 <- gasMatNz[fdr20.idx,]
bc20.fdr20.sqrt <- bcSubSamplePar(sqrt(fdr20), 100, 20, 0.6)
bc25.fdr20.sqrt <- bcSubSamplePar(sqrt(fdr20), 100, 25, 0.6)
bc30.fdr20.sqrt <- bcSubSamplePar(sqrt(fdr20), 100, 30, 0.6)

fdr30p90Z <- t(scale(t(fdr30p90)))
bc30.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 30, 0.6)
bc25.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 25, 0.6)
bc20.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 20, 0.6)

save(bc30.fdr30p90Z, file = "bc30.fdr30p90Z.RData")

fdr20Z <- t(scale(t(fdr20)))
bc20.fdr20Z <- bcSubSamplePar(fdr20Z, 100, 20, 0.6)
bc25.fdr20Z <- bcSubSamplePar(fdr20Z, 100, 25, 0.6)
bc30.fdr20Z <- bcSubSamplePar(fdr20Z, 100, 30, 0.6)
save(bc20.fdr20Z,
     bc25.fdr20Z,
     bc30.fdr20Z, file = "bcFdr20Z.RData")


save(bc30.fdr30p90.sqrt,
     bc25.fdr30p90.sqrt,
     bc20.fdr30p90.sqrt,
     bc20.fdr20.sqrt,
     bc25.fdr20.sqrt,
     bc30.fdr20.sqrt,
     bc30.fdr30p90Z,
     bc25.fdr30p90Z,
     bc20.fdr30p90Z,
     bc20.fdr20Z,
     bc25.fdr20Z,
     bc30.fdr20Z, file = "bcGas.RData")

# making figures...
load("../bcSol/bcGas.RData")

datTmp <- t(scale(t(fdr30p90)))
pssTmp <- postSubSample.percent(bc30.fdr30p90.sqrt, 0.9, 0.5)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

datTmp <- sqrt(sqrt(fdr30p90))

datTmp <- t(scale(t(fdr30p90)))
pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.9, 0.5)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

# this bicluster looks good... but lots of zeroes... get rid of zeroes and find next BC?
load("../bcSol/bcFdr20Z.RData")
datTmp <- sqrt(fdr20)

datTmp <- (fdr20Z)
pssTmp <- postSubSample.percent(bc30.fdr20Z , 0.9, 0.5)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
write(rownames(fdr20Z[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr20Z_FG.txt")
write(rownames(fdr20Z), file = "../results/fdr20Z_BG.txt")

clustIdxTmp <- get1DIdx.ext(fdr20Z,pssTmp$rowIdx, pssTmp$colIdx)
fdr20.minus1 <- fdr20
fdr20.minus1[clustIdxTmp] <- sample(fdr20[-clustIdxTmp], length(clustIdxTmp))
fdr20.minus1Z <- t(scale(t(fdr20.minus1)))

bc30.fdr20Z.2 <- bcSubSamplePar(fdr20.minus1Z, 100, 30, 0.6)
save(bc30.fdr20Z.2, file = "bcFdr20Z.2.RData")

load("../bcSol/bcFdr20Z.2.RData")
datTmp <- fdr20Z
datTmp <- sqrt(fdr20)
pssTmp <- postSubSample.percent(bc30.fdr20Z.2 , 0.65, 0.8)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
write(rownames(datTmp[pssTmp$rowIdx,]), file = "../results/fdr20Z_2_FG.txt")

pssTmp <- postSubSample.percent(bc20.fdr20Z , 0.75, 0.9)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

# look at 2nd w/ 20 conditions
# 25 seems to be a better condition reg
pssTmp <- postSubSample.percent(bc25.fdr20Z , 0.75, 0.9)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggsave("../img/fdr20Z_75_90.pdf", width = 21.6, height = 9.91)

write(rownames(datTmp[pssTmp$rowIdx,]), file = "../results/fdr20Z_25_FG.txt") # this looks pretty terrible

clustIdxTmp <- get1DIdx.ext(fdr20Z,pssTmp$rowIdx, pssTmp$colIdx)
fdr20.minus1 <- fdr20
fdr20.minus1[clustIdxTmp] <- sample(fdr20[-clustIdxTmp], length(clustIdxTmp))
fdr20.minus1Z <- t(scale(t(fdr20.minus1)))

bc25.fdr20Z.25.2 <- bcSubSamplePar(fdr20.minus1Z, 100, 25, 0.6)
bc25.fdr20Z.20.2 <- bcSubSamplePar(fdr20.minus1Z, 100, 20, 0.6)
save(bc25.fdr20Z.25.2,
     bc25.fdr20Z.20.2, file = "bc25.2.RData")


load("../bcSol/bc25.2.RData")

pssTmp <- postSubSample.percent(bc25.fdr20Z.25.2, 0.75, 0.75)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])



# looking at FDR and no zeroes
load("../bcSol/bc30.fdr30p90Z.RData")

datTmp <- (fdr30p90)
datTmp <- t(scale(t(fdr30p90)))
datTmp <- sqrt(sqrt(fdr30p90))
pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.9, 0.5)

pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.90, 0.90)
write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr30p90Z_90_90_FG.txt")

pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.75, 0.90)
write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr30p90Z_75_90_FG.txt")

write.table(datTmp[pssTmp$rowIdx, pssTmp$colIdx], file = "../results/fdr30p90Z_75_90_reactome_raw.txt",
            sep = "\t", col.names = F, quote = F)
write.table(datTmp[pssTmp$rowIdx, pssTmp$colIdx], file = "../results/fdr30p90Z_75_90_reactome_z.txt",
            sep = "\t", col.names = F, quote = F)


write(colnames(datTmp[, pssTmp$colIdx]), file = "../resultsForBen/conditionsInBC.txt")

# FDR 30, 90% of conditions don't have zeroes
# Fig. to show Prof. Huang
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggsave("../img/fdr30p90Z_75_90.pdf", width = 21.6, height = 9.91)

# 75 90 looks good... lets find 2nd BC
pssTmp <- postSubSample.percent(bc30.fdr30p90Z , 0.75, 0.90)
fdr30p90.minus1 <- maskData(fdr30p90, pssTmp)
fdr30p90.minus1Z <- t(scale(t(fdr30p90.minus1)))
bc25.fdr30p90.minus1Z <- bcSubSamplePar(fdr30p90.minus1Z, 100, 25, 0.6)
bc30.fdr30p90.minus1Z <- bcSubSamplePar(fdr30p90.minus1Z, 100, 30, 0.6)
save(bc25.fdr30p90.minus1Z,
     bc30.fdr30p90.minus1Z, file = "bc2.fdr30p90.RData")

load("../bcSol/bc2.fdr30p90.RData")
pssTmp <- postSubSample.percent(bc25.fdr30p90.minus1Z , 0.75, 0.75)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

pssTmp <- postSubSample.percent(bc30.fdr30p90.minus1Z , 0.75, 0.85)
write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr30p90Z.2_FG.txt")
tmpExpression <- datTmp[pssTmp$rowIdx, pssTmp$colIdx]
corTmp <- (cor(t(tmpExpression)))
hist(corTmp[upper.tri(corTmp)])
hist(cor(t(tmpExpression)))
hist(cor(t(flyWorm)))


ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr30p90Z_FG.txt")
write(rownames(datTmp), file = "../results/fdr30p90Z_BG.txt")

load("../bcSol/bcFdr30p90.RData")
datTmp <- t(scale(t(fdr30p90)))
datTmp <- sqrt(fdr30p90)
datTmp <- (fdr30p90)
datTmp <- sqrt(sqrt(fdr30p90))
pssTmp <- postSubSample.percent(bc30.fdr30p90.sqrt , 0.9, 0.5)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotExpression(datTmp[sample(1:nrow(datTmp), 180), pssTmp$colIdx])
write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../results/fdr30p90sqrt_FG.txt")


# finding smaller biclusters
bc5.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 5, 0.6)
bc8.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 8, 0.6)
bc10.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 100, 10, 0.6)
