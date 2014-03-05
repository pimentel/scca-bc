library(gplots)

flyWorm <- read.table("../data/flyWorm.txt", header = TRUE)
flyWorm <- as.matrix(flyWorm)

fwBC <- bcSubSamplePar(flyWorm, 100, 27, 0.6)
save(fwBC, file = "fwBC.RData")

fwBC30 <- bcSubSamplePar(flyWorm, 100, 30, 0.6)


fwBC33 <- bcSubSamplePar(flyWorm, 100, 33, 0.6)

fwBC35 <- bcSubSamplePar(flyWorm, 100, 35, 0.6)

fwBC27.all <- biclusteringPar(flyWorm, 100, 27)

fwBC30.all <- biclusteringPar(flyWorm, 100, 30)

fwBC33.all <- biclusteringPar(flyWorm, 100, 33)

fwBC35.all <- biclusteringPar(flyWorm, 100, 35)

save(fwBC30, fwBC33, fwBC35, fwBC27.all, fwBC30.all, fwBC33.all, 
     fwBC35.all, file = "fwBCExtra.RData")

load("fwBCExtra.RData", verbose = T)
tmpAllNames <- rownames(flyWorm)
rownames(flyWorm) <- flybaseNames

pss.fwBC30 <- postSubSample.percent(fwBC30, 0.95, 0.5)

pss.fwBC30 <- postSubSample.percent(fwBC30, 0.9, 0.9)
write((flybaseNames[pss.fwBC30$rowIdx]), file = "../results/flybaseNames_30_90_90.FG.txt") # this is most interesting
write((wormNames[pss.fwBC30$rowIdx]), file = "../results/wormNames_30_90_90.FG.txt") # this is most interesting


write.table(flyWorm[pss.fwBC30$rowIdx, pss.fwBC30$colIdx], 
            file = "../results/flybaseNames_30_90_90_reactome.txt", 
            sep = "\t", col.names = F, quote = F)

rownames(flyWorm) <- wormNames
write.table(flyWorm[pss.fwBC30$rowIdx, pss.fwBC30$colIdx], 
            file = "../results/wormNames_30_90_90_reactome.txt", 
            sep = "\t", col.names = F, quote = F)

nrow(flyWorm[pss.fwBC30$rowIdx, pss.fwBC30$colIdx])

pss.fwBC30 <- postSubSample.percent(fwBC30, 0.95, 0.95)
write((flybaseNames[pss.fwBC30$rowIdx]), file = "../results/flybaseNames_30_95_95.FG.txt") # this is most interesting
write((flybaseNames[pss.fwBC30$rowIdx]), file = "../results/flybaseNames_30_95_95.FG.txt") # this is most interesting
colnames(flyWorm)[pss.fwBC30$colIdx]
pss.fwBC30$colIdx

maskData <- function(geneDf, clusterIdx)
{
    # getting next BC
    get1DIdx.ext <- function(geneDf, rows, cols) 
    {
        allPairs <- expand.grid(rows, cols)
        nrow(geneDf) * (allPairs[, 2] - 1) + allPairs[, 1]
    }
    idxToMask <- get1DIdx.ext(geneDf, clusterIdx$rowIdx, clusterIdx$colIdx)
    geneDf[idxToMask] <- sample(geneDf[-idxToMask], length(idxToMask))
    return(geneDf)
}

# 0.9, 0.9 looks good... lets get next bicluster
pss.fwBC30 <- postSubSample.percent(fwBC30, 0.9, 0.9)
flyWorm.minus1 <- maskData(flyWorm, pss.fwBC30)
fw30.2 <- bcSubSamplePar(flyWorm.minus1, 100, 30, 0.6)
save(fw30.2, file = "fw30.2.RData")

fw25.2 <- bcSubSamplePar(flyWorm.minus1, 100, 25, 0.6)
save(fw25.2, file = "fw25.2.RData")

load("../bcSol/fw30.2.RData")
pssTmp <- postSubSample.percent(fw30.2, 0.75, 0.9)
write((flybaseNames[pssTmp$rowIdx]), file = "../results/flybaseNames_30_75_90_2.FG.txt") # this is most interesting
write((wormNames[pssTmp$rowIdx]), file = "../results/wormNames_30_75_90_2.FG.txt") # this is most interesting
pssTmp <- postSubSample.percent(fw30.2, 0.90, 0.9)
write((flybaseNames[pssTmp$rowIdx]), file = "../results/flybaseNames_30_90_90_2.FG.txt") # this is most interesting
write((wormNames[pssTmp$rowIdx]), file = "../results/wormNames_30_90_90_2.FG.txt") # this is most interesting
ggPlotExpression(flyWorm[pssTmp$rowIdx, pssTmp$colIdx])

ggPlotParSolution(fw30.2)





set.seed(42)
ranRows <- sample(setdiff(1:nrow(flyWorm), pss.fwBC30$rowIdx), 200)
fwBC30.trim <- lapply(fwBC30, function(x) {
                      x$ab <- x$ab[c(pss.fwBC30$rowIdx, ranRows)]
                      x
     })

plotParSolution(fwBC30.trim, "../img/flyWorm")


A <- abs(sapply(fwBC30.trim, function(x) x$ab))

ggPlotParSolution(fwBC30.trim, "flyWorm", colNames = colnames(flyWorm))


plotParSolution(fwBC30)

levelplot(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx))], 
          scales=list(x=list(rot=90)))
dev.print(pdf, file = "flyWorm30.pdf")

colnames(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx))])

ggPlotExpression(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx))])
ggPlotExpression(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx))])

# TODO: fix this figure
# Fly/worm Fig. to show Prof. Huang
ggPlotExpression(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx))])
ggsave("../img/flyWorm_90_90.pdf", width = 21.6, height = 9.91 )
# Saving 21.6 x 9.91 in image
 

levelplot(flyWorm[pss.fwBC30$rowIdx, sort(c(pss.fwBC30$colIdx, 32))], 
          scales=list(x=list(rot=90)))
dev.print(pdf, file = "flyWorm30PlusLate.pdf")



pss.fwBC33 <- postSubSample.percent(fwBC33, 0.95, 0.5)
colnames(flyWorm)[pss.fwBC33$colIdx]
pss.fwBC33$colIdx
levelplot(flyWorm[pss.fwBC33$rowIdx, sort(c(pss.fwBC33$colIdx))], 
          scales=list(x=list(rot=90)))

pss.fwBC27.all <- postSubSample.percent(fwBC27.all, 0.95, 0.5)
colnames(flyWorm)[pss.fwBC27.all$colIdx]
pss.fwBC27.all$colIdx
levelplot(flyWorm[pss.fwBC27.all$rowIdx, sort(c(pss.fwBC27.all$colIdx))], 
          scales=list(x=list(rot=90)))

pss.fwBC30.all <- postSubSample.percent(fwBC30.all, 0.85, 0.5)
colnames(flyWorm)[pss.fwBC30.all$colIdx]
pss.fwBC30.all$colIdx
levelplot(flyWorm[pss.fwBC30.all$rowIdx, sort(c(pss.fwBC27.all$colIdx))], 
          scales=list(x=list(rot=90)))

post.fwBC <- postSubSample(fwBC)

post.fwBC.hclust <- post.hclust(post.fwBC, 3, 2)
post.fwBC.hclust

post.iid.subSample.per <- postSubSample.percent(fwBC, 0.95, .5)
post.iid.subSample.per

hi <- sapply(fwBC, function(x) x$d)
test2 <- postSubSample.percent2(fwBC, 0.75, .3)


order.fw <- flyWorm[post.iid.subSample.per$rowIdx, post.fwBC.hclust[[1]]$colIdx]
levelplot(order.fw, scales=list(x=list(rot=90)))

flybaseNames <- sapply(strsplit(rownames(flyWorm), "/"), function(x) x[1])
wormNames <- sapply(strsplit(rownames(flyWorm), "/"), function(x) x[2])
write(wormNames, "../results/wormNames_bg.txt")
write(flybaseNames, "../results/flybaseNames_BG.txt")
write(flybaseNames[post.iid.subSample.per$rowIdx], "../results/flybaseNames_FG.txt")

write(flybaseNames[pss.fwBC30$rowIdx], "../results/flybaseNames_30_FG.txt")

write(flybaseNames[pss.fwBC30$rowIdx], "../results/flybaseNames_30_85p_FG.txt")

dev.print(pdf, "../img/flyWorm.pdf")

levelplot(order.fw, col.regions = redgreen(75),
          scales=list(x=list(rot=90)))

geneNames <- simplify2array(strsplit(rownames(flyWorm)[post.iid.subSample.per$rowIdx], "/"))[2,]
levelplot(flyWorm[post.fwBC.hclust[[1]]$rowIdx, post.fwBC.hclust[[1]]$colIdx],
          scales=list(x=list(rot=90)), )

levelplot(flyWorm[sample(nrow(flyWorm), 148), (sample(ncol(flyWorm), 27))],
          scales=list(x=list(rot=90)))

dev.print(pdf, "../img/flyWormRandom.pdf")

levelplot(flyWorm[sample(nrow(flyWorm), 148), post.fwBC.hclust[[1]]$colIdx],
          scales=list(x=list(rot=90)))



load("fwBC.RData")
plotParSolution(fwBC)


post.fwBC.kmeans <- post.hclust(post.fwBC, 3, 2)
post.fwBC.kmeans


# go terms for worm
sexDiff<-c("C04H5.6","F56G4.4","W03F9.10","B0286.4","W09C5.2","R07E5.3","K07A1.12","C14B9.4","F41E6.4","C01G8.9","T10F2.4","Y116A8C.32","Y41D4B.19","R05D3.4","F56D2.6","W07E6.4","W04A8.7","ZK507.6","R53.6","Y110A7A.8","F26E4.10","C43E11.10","T13H5.4","T05G5.3","C06A8.2","C08B11.3","JC8.6","Y92H12A.1","Y57A10A.19","C36B1.5","F55F8.4","T23G7.1","W02D3.9","F59E10.2","F10B5.6","F12F6.3","Y111B2A.22","C55A6.9","T12A2.7","W02A11.4","R08D7.1","B0035.11","Y54E5B.3","T11G6.8","ZK616.4","F58A4.4","C15C6.4","F25B3.6","Y113G7B.23","F19F10.9","K08E4.1","F09G2.4","T08A11.2","Y106G6E.5","M04B2.1","C50C3.6")
herm<-c("C04H5.6","F56G4.4","W03F9.10","B0286.4","W09C5.2","R07E5.3","K07A1.12","C14B9.4","F41E6.4","C01G8.9","T10F2.4","Y116A8C.32","Y41D4B.19","R05D3.4","F56D2.6","W07E6.4","W04A8.7","ZK507.6","R53.6","Y110A7A.8","C43E11.10","T13H5.4","T05G5.3","C06A8.2","JC8.6","C08B11.3","Y57A10A.19","C36B1.5","F55F8.4","T23G7.1","W02D3.9","F10B5.6","F12F6.3","Y111B2A.22","C55A6.9","W02A11.4","T12A2.7","R08D7.1","B0035.11","Y54E5B.3","T11G6.8","ZK616.4","F58A4.4","C15C6.4","F25B3.6","Y113G7B.23","F19F10.9","K08E4.1","F09G2.4","T08A11.2","M04B2.1","C50C3.6")
gen<-c("C04H5.6","F56G4.4","W03F9.10","B0286.4","W09C5.2","R07E5.3","K07A1.12","C14B9.4","F41E6.4","C01G8.9","T10F2.4","Y116A8C.32","Y41D4B.19","R05D3.4","F56D2.6","W07E6.4","W04A8.7","ZK507.6","R53.6","Y110A7A.8","C43E11.10","T13H5.4","T05G5.3","C06A8.2","JC8.6","C08B11.3","Y57A10A.19","C36B1.5","F55F8.4","T23G7.1","W02D3.9","F10B5.6","F12F6.3","Y111B2A.22","C55A6.9","W02A11.4","T12A2.7","R08D7.1","B0035.11","Y54E5B.3","T11G6.8","ZK616.4","F58A4.4","C15C6.4","F25B3.6","Y113G7B.23","F19F10.9","K08E4.1","F09G2.4","T08A11.2","M04B2.1","C50C3.6")


clust2<-c("F28C6.3","C04H5.6","C50F2.3","E01A2.4","F49D11.1","D1081.8","Y54E5B.3","C36B1.3","Y102A5C.18","C36B1.5","F55F8.4","T05H4.14","F33A8.1","K01G5.6","R10E4.4","W04A8.7","Y54E10BR.6","K07A1.12","F12F6.3","C01G8.9","Y110A7A.8")
clust2.2 <- c("F28C6.3","C50F2.3","E01A2.4","F49D11.1","D1081.8","Y54E5B.3","C36B1.3","C36B1.5","F55F8.4","T05H4.14","F33A8.1","R10E4.4","W04A8.7","Y54E10BR.6","C01G8.9","Y110A7A.8")


# out of curiosity, look at smaller BC with 13 conditions
fwBC13 <- bcSubSamplePar(flyWorm, 100, 13, 0.6)
save(fwBC13, file = "fwBC13.RData")

fwBC10 <- bcSubSamplePar(flyWorm, 100, 10, 0.6)
pssTmp <- postSubSample.percent(fwBC10, 0.9, 0.9)


load("../bcSol/fwBC13.RData")

pssTmp <- postSubSample.percent(fwBC13, 0.9, 0.9)
ggPlotExpression(flyWorm[pssTmp$rowIdx, pssTmp$colIdx])

