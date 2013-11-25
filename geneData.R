geneData <- read.csv("../data/Supplementary_Data_File_10_gene_expression_BPKM.csv",
                     header = T)

geneDataAdult <- geneData[, -grep("L3", colnames(geneData))]
geneDataAdult <- geneData[, order(colnames(geneDataAdult))]
gas <- geneDataAdult[, order(colnames(geneDataAdult))]
gas.mat <- as.matrix(gas)

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

save(gasBC.raw, file = "gasBC.raw.RData")

load("gasBC.raw.RData")

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
sum(which.20.fdr)

which.20.fdr[which(!(which.20.fdr %in% rownames(logGA)))]
validIdx <- which.20.fdr[which((which.20.fdr %in% rownames(gas.mat)))]
