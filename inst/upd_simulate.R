cm <- constMean.norm(1)

load("~/Dropbox/biclustering/R/constMean.sol.RData", verbose = T)

data <- constMean.sol$data[[1]]$mat

# res <- bcSubSampleSerial(data, 50, 17)
res <- biclusteringSerial(data, 50, 15)

res <- bcSubSampleSerial(data, 50, 15)
x <- ggPlotSSSolution(res)

a <- postSubSample.pca(res)

load("~/Dropbox/biclustering/R/sim50Mvn.sol.RData", verbose = T)

data <- sim50Mvn.sol$data[[1]]$mat
res <- biclusteringSerial(data, 50, 30)

x <- ggPlotSSSolution(res)

debugonce(postSubSample.pca)
a <- postSubSample.pca(res)

debugonce(postSubSample.ranks)
r <- postSubSample.ranks(res)

boxplot(r)
ggplot(data.frame(x = 1:1500, y = apply(r, 1, median)), aes(x, y)) + geom_point()

jaccard_idx_matrix(list(r), list(list(rowIdx = 1:300, colIdx = 1:30)))

jaccard_idx_matrix(list(a), list(list(rowIdx = 1:300, colIdx = 1:30)))

jaccard_idx_matrix(sim50Mvn.sol$sol[[1]]$clusters, list(list(rowIdx = 1:300, colIdx = 1:30)))

prev_sol <- postSubSample.pca(sim50Mvn.sol$sol[[1]]$allBcSol[[1]])
jaccard_idx_matrix(list(prev_sol), list(list(rowIdx = 1:300, colIdx = 1:30)))

ggplot(r$ab_features, aes(PC1, PC2, group = in_clust, colour = in_clust)) + geom_point()

ggplot(r$d_features, aes(PC1, PC2, group = in_clust, colour = in_clust)) + geom_point()

ggplot(a$abDat, aes(-PC1 * sdev[1], PC2 * sdev[2], colour = cluster)) + geom_point()

ggplot(a$abDat, aes(PC1, PC2)) + geom_point()

ggplot(data.frame(a$dDat$rotation), aes(PC1, PC2)) + geom_point()

x <- ggPlotSSSolution(res)
y <- ggPlotSSSolution(sim50Mvn.sol$sol[[1]]$allBcSol[[1]])

grid.arrange(x$condPlot, y$condPlot + ylim(0, 50))
