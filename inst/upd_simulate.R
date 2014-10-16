
cm <- constMean.norm(1)

load("~/Dropbox/biclustering/R/constMean.sol.RData", verbose = T)

data <- constMean.sol$data[[1]]$mat

res <- bcSubSampleSerial(data, 50, 17)
x <- ggPlotSSSolution(res)

