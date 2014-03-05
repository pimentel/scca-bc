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

# 2nd bicluster...not particularly interesting... mostly adult stages
load("../bcSol/bc2.fdr30p90.RData")
pssTmp <- postSubSample.percent(bc25.fdr30p90.minus1Z , 0.75, 0.75)
pssTmp <- postSubSample.top(bc25.fdr30p90.minus1Z , 0.75, 0.75, 25)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggsave("../resultsBen2/2nd_bicluster.pdf")

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
# UPDATE: these all look like crap...
bc5.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 200, 5, 0.6)
bc8.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 200, 8, 0.6, lam.lwr = 5)
bc10.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 200, 10, 0.6, lam.lwr = 7)
bc15.fdr30p90Z <- bcSubSamplePar(fdr30p90Z, 200, 15, 0.6, lam.lwr = 7)
save(bc5.fdr30p90Z, 
     bc8.fdr30p90Z,
     bc10.fdr30p90Z, 
     bc15.fdr30p90Z, file = "bcSmall.fdr30.RData")
save(bc5.fdr30p90Z, file = "bc5.fdr30p90Z.RData")


load("../bcSol/bc5.fdr30p90Z.RData")
postSubSample.percent(bc5.fdr30p90Z, 0.9, 0.5)
pssTmp <- postSubSample.percent2(bc5.fdr30p90Z, 0.85, 0.5, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

load("../bcSol/bcSmall.fdr30.RData")
pssTmp <- postSubSample.percent2(bc5.fdr30p90Z, 0.85, 0.6, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

tmpCor <- cor(t(datTmp[pssTmp$rowIdx, pssTmp$colIdx]))
hist(tmpCor[upper.tri(tmpCor)], xlim = c(-1, 1))

pssTmp <- postSubSample.percent(bc8.fdr30p90Z, 0.85, 0.6)
pssTmp <- postSubSample.percent2(bc8.fdr30p90Z, 0.85, 0.6, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])


pssTmp <- postSubSample.percent(bc10.fdr30p90Z, 0.85, 0.6)
pssTmp <- postSubSample.percent2(bc10.fdr30p90Z, 0.85, 0.6, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

pssTmp <- postSubSample.percent(bc15.fdr30p90Z, 0.85, 0.6)
pssTmp <- postSubSample.percent2(bc15.fdr30p90Z, 0.85, 0.6, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

# look at samples that Ben is interested in 
benSamples <- (c("Zn_4.5mM_AdMixed_4days","Cd_0.05M_Cdcl2_AdMixed_4days","Cd_0.1M_Cdcl2_AdMixed_4days","Cu_15mM_AdMixed_4days","Cold1_AdMixed_4days","Cold2_AdMixed_4days","HeatShock1hr36C_WetHeat30min25C_AdMixed_4days","Caffeine_2.5mg_per_ml_AdMixed_4days","Caffeine_PreLethal_AdMixed_4days","Paraquat_5mM_AdMixed_4days","Paraquat_10mM_AdMixed_4days"))
benIdx <- trimPercent(gas.mat[,benSamples], 
                      idxAtFDR(geneDE.fdr, rownames(gas.mat), 0.3), 0.9)
benMat <- gas.mat[benIdx, benSamples]
benMatZ <- t(scale(t(benMat)))

bc5.ben <- bcSubSamplePar(benMatZ, 200, 5, 0.6, lam.lwr = 4)
bc7.ben <- bcSubSamplePar(benMatZ, 200, 7, 0.6, lam.lwr = 5)
bc9.ben <- bcSubSamplePar(benMatZ, 200, 9, 0.6, lam.lwr = 5)
bc11.ben <- bcSubSamplePar(benMatZ, 200, 11, 0.6, lam.lwr = 5)
save(bc5.ben, 
     bc7.ben,     
     bc9.ben,     
     bc11.ben,file = "bc5.ben.RData")

load("../bcSol/bc5.ben.RData")

# this is kind of interesting... up and down regulated genes
datTmp <- benMatZ
datTmp <- benMat
pssTmp <- postSubSample.percent(bc5.ben, 0.9, 0.5)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggsave("../resultsBen2/perturb_90_50.pdf")
write(rownames(datTmp[pssTmp$rowIdx, pssTmp$colIdx]), file = "../resultsBen2/perturb5Z_75_90_FG.txt")
write(rownames(datTmp), file = "../resultsBen2/perturb_BG.txt")
write.table(datTmp[pssTmp$rowIdx, pssTmp$colIdx], file = "../resultsBen2/perturbZ_90_50_reactome_z.txt",
            sep = "\t", col.names = F, quote = F)


postSubSample.percent2(bc5.ben, 0.9, 0.5, 0.9, 0.5)
postSubSample.top(bc5.ben, 0.9, 0.5, 5)

pssTmp <- postSubSample.percent(bc7.ben, 0.9, 0.5)
pssTmp <- postSubSample.percent(bc9.ben, 0.85, 0.9)
pssTmp <- postSubSample.top(bc9.ben, 0.85, 0.9, 9)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc9.ben)

pssTmp <- postSubSample.percent(bc11.ben, 0.95, 0.9)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
tmpCor <- cor(t(datTmp[pssTmp$rowIdx, pssTmp$colIdx]))
hist(tmpCor[upper.tri(tmpCor)])


# next... try to remove all the "normal" samples
normNames <- grep("^Ad", colnames(gas.mat))
notNorm <- setdiff(1:ncol(gas.mat), normNames)
notNormIdx <- trimPercent(gas.mat[,notNorm], 
                          idxAtFDR(geneDE.fdr, rownames(gas.mat), 0.3), 0.9)
nnMat <- gas.mat[notNormIdx, notNorm]
nnZ <- t(scale(t(nnMat)))

bc5.nn <- bcSubSamplePar(nnZ, 200, 5, 0.6, lam.lwr = 4)
bc7.nn <- bcSubSamplePar(nnZ, 200, 7, 0.6, lam.lwr = 5)
bc9.nn <- bcSubSamplePar(nnZ, 200, 9, 0.6, lam.lwr = 5)
bc11.nn <- bcSubSamplePar(nnZ, 200, 11, 0.6, lam.lwr = 5)
bc15.nn <- bcSubSamplePar(nnZ, 200, 11, 0.6, lam.lwr = 7)
save(bc5.nn, 
     bc7.nn, 
     bc9.nn, 
     bc11.nn,
     bc15.nn,file = "bc.nn.RData")

load("../bcSol/bc.nn.RData")

datTmp <- nnZ

pssTmp <- postSubSample.percent(bc5.nn, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc5.nn, 0.75, 0.75, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

pssTmp <- postSubSample.percent(bc7.nn, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc7.nn, 0.75, 0.75, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

pssTmp <- postSubSample.percent(bc9.nn, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc9.nn, 0.75, 0.75, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])


pssTmp <- postSubSample.percent(bc11.nn, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc11.nn, 0.75, 0.75, 0.75, 0.50)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

pssTmp <- postSubSample.percent(bc15.nn, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc15.nn, 0.75, 0.75, 0.5, 0.50)
ggPlotParSolution(bc15.nn)
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])

par(mfrow = c(2, 1))
tmpCor <- cor(t(datTmp[pssTmp$rowIdx, pssTmp$colIdx]))
hist(tmpCor[upper.tri(tmpCor)])

# remove embryo stuff too
nAE <- grep("^Ad|^em", colnames(gas.mat))
nAECols <- setdiff(1:ncol(gas.mat), nAE)
nAERows <- trimPercent(gas.mat[,nAECols],
                       idxAtFDR(geneDE.fdr, rownames(gas.mat), 0.3),
                       0.945)
nAEMat <- gas.mat[nAERows, nAECols]
nAEMatZ <- t(scale(t(nAEMat)))

datTmp <- nAEMatZ
datTmp <- nAEMat

bc5nAE <- bcSubSamplePar(nAEMatZ, 200, 5, 0.6, lam.lwr = 4)
bc7nAE <- bcSubSamplePar(nAEMatZ, 200, 7, 0.6, lam.lwr = 5)
bc9nAE <- bcSubSamplePar(nAEMatZ, 200, 9, 0.6, lam.lwr = 5)
bc11nAE <- bcSubSamplePar(nAEMatZ, 200, 11, 0.6, lam.lwr = 5)
bc15nAE <- bcSubSamplePar(nAEMatZ, 200, 15, 0.6, lam.lwr = 12)
save(bc5nAE,
     bc7nAE,
     bc9nAE,
     bc11nAE,
     bc15nAE, file = "bc.nAE.RData")

load("../bcSol/bc.nAE.RData")

pssTmp <- postSubSample.percent(bc5nAE, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc5nAE, 0.75, 0.75, 0.60, 0.50)
pssTmp <- postSubSample.top(bc5nAE, 0.75, 0.75, 5)
pssTmp
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc5nAE)

pssTmp <- postSubSample.percent(bc7nAE, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc7nAE, 0.75, 0.75, 0.60, 0.50)
pssTmp
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc5nAE)


pssTmp <- postSubSample.percent(bc9nAE, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc9nAE, 0.75, 0.75, 0.60, 0.50)
pssTmp
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc9nAE)


pssTmp <- postSubSample.percent(bc11nAE, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc11nAE, 0.75, 0.75, 0.60, 0.50)
pssTmp
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc11nAE)

pssTmp <- postSubSample.percent(bc15nAE, 0.75, 0.75)
pssTmp <- postSubSample.percent2(bc15nAE, 0.80, 0.75, 0.50, 0.50)
pssTmp <- postSubSample.sort(bc15nAE, 0.75, 0.75)
pssTmp
ggPlotExpression(datTmp[pssTmp$rowIdx, pssTmp$colIdx])
ggPlotParSolution(bc15nAE)



