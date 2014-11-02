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
colnames(data) <- sapply(1:ncol(data), function(x) paste(sample(letters, 3), collapse = ""))
params <- sccab_params(30, n_samp = 50, ab_lam = c(1, 10, 20))
res <- sccab(data, params)
sr <- sccab_result(res, params)

params <- sccab_params(30, n_samp = 50, ab_lam = c(1, 10, 20), prop = 0.6, parallel = F)
res1 <- sccab_subsample(data, params)
hi <- sccab_result(res1, FALSE,params)
x <- ggPlotSSSolution(res1)

res2 <- sccab(data, 30, 50, clust_opt = list(lamx = c(1, 2, 3)))

x <- ggPlotSSSolution(res)

debugonce(postSubSample.pca)
a <- postSubSample.pca(res)
a1 <- postSubSample.pca(res1)
a2 <- postSubSample.pca(res2)

debugonce(postSubSample.ranks)
r <- postSubSample.ranks(res)
r1 <- postSubSample.ranks(res1)
r2 <- postSubSample.ranks(res2)

h <- post_hclust(res)
h2 <- post_hclust(res2)

boxplot(r)
ggplot(data.frame(x = 1:1500, y = apply(r, 1, median)), aes(x, y)) + geom_point()

jaccard_idx_matrix(list(r), list(list(rowIdx = 1:300, colIdx = 1:30)))
jaccard_idx_matrix(list(r1), list(list(rowIdx = 1:300, colIdx = 1:30)))
jaccard_idx_matrix(list(r2), list(list(rowIdx = 1:300, colIdx = 1:30)))

jaccard_idx_matrix(list(h), list(list(rowIdx = 1:300, colIdx = 1:30)))
jaccard_idx_matrix(list(h2), list(list(rowIdx = 1:300, colIdx = 1:30)))

jaccard_idx_matrix(list(a), list(list(rowIdx = 1:300, colIdx = 1:30)))
jaccard_idx_matrix(list(a1), list(list(rowIdx = 1:300, colIdx = 1:30)))
jaccard_idx_matrix(list(a2), list(list(rowIdx = 1:300, colIdx = 1:30)))

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
