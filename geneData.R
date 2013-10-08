geneData <- read.csv("../data/Supplementary_Data_File_10_gene_expression_BPKM.csv",
                     header = T)

geneDataAdult <- geneData[, -grep("L3", colnames(geneData))]
geneDataAdult <- geneData[, order(colnames(geneDataAdult))]
gas <- geneDataAdult[, order(colnames(geneDataAdult))]
gas.mat <- as.matrix(gas)

d <- diag(runif(ncol(gas)))
colnames(d) <- colnames(gas.mat)



