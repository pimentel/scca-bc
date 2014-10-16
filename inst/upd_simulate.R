
cm <- constMean.norm(1)

load("~/Dropbox/biclustering/R/constMean.sol.RData", verbose = T)

data <- constMean.sol$data[[1]]$mat

res <- bcSubSampleSerial(data, 50, 17)
x <- ggPlotSSSolution(res)


load("~/Dropbox/biclustering/R/sim50Mvn.sol.RData", verbose = T)

data <- sim50Mvn.sol$data[[1]]$mat
res <- biclusteringSerial(data, 50, 30)
x <- ggPlotSSSolution(res)
y <- ggPlotSSSolution(sim50Mvn.sol$sol[[1]]$allBcSol[[1]])

grid.arrange(x$condPlot, y$condPlot + ylim(0, 50))
