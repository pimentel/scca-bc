library(vioplot)

# First simulation. Basic constant mean bicluster

constMean.norm <- function(nSamples)
{
    set.seed(42)

    blocks <- vector(nSamples, mode = "list")
    solz <- vector(nSamples, mode = "list")

    for (it in 1:nSamples )
    {
        cat("***********************************\n")
        cat("Running simulation: ", it, "\n")
        cat("***********************************\n")
        curBicluster <- genBlock.iid.N(30, 15, 4, 1)
        curSim <- generateNormal(nrow = 300, ncol = 40, clusterOptions = 
                               list(list(x.start = 1, x.end = 30, y.start = 1, y.end = 15,
                                         data = curBicluster)))
        ps.curSim <- permuteSim(curSim)
        blocks[[it]] <- ps.curSim
        curSol <- bcMultipleClusters(ps.curSim$permutedMat, 15, 1, 100)
        curSol.reorder <- reorderBCSol(curSol, ps.curSim)
        solz[[it]] <- curSol.reorder
    }


    list(data = blocks, sol = solz)
}

computePerf <- function(truth, solList, metric)
{
    sapply(solList, function(sol)
           {
               metric(truth, sol$clusters)
           })
}

debug(constMean.norm)
undebug(constMean.norm)

constMean.sol <- constMean.norm(50)
save(constMean.sol, file = "constMean.sol.RData")
constMean.truth <- list(list(rowIdx = 1:30, colIdx = 1:15))

load("constMean.sol.RData")

constMean.perf <- computePerf(constMean.truth, constMean.sol$sol, amrs.hp)


computePerf(constMean.truth, constMean.sol$sol, amrs)
mean(computePerf(constMean.truth, constMean.sol$sol, amrs.hp))
median(computePerf(constMean.truth, constMean.sol$sol, amrs.hp))
sd(computePerf(constMean.truth, constMean.sol$sol, amrs.hp))

# performance sim 1
plot(c(1,1), c(0, 1), xlim = c(1, 1), ylim = c(0, 1))
vioplot(constMean.perf,
        col = "lightblue", add = T
        )
boxplot(constMean.perf, ylim = c(0, 1))

# look at a good one and bad one
plotParSolution(constMean.sol$sol[[1]]$allBcSol[[1]], "constMeanGood")
ggPlotParSolution(constMean.sol$sol[[1]]$allBcSol[[1]], "constMeanGood")
plotParSolution(constMean.sol$sol[[20]]$allBcSol[[1]], "constMeanBad")
ggPlotParSolution(constMean.sol$sol[[20]]$allBcSol[[1]], "constMeanBad")


simSol <- constMean.norm(50)
constMean.sol <- simSol
save(constMean.sol, file = "constMean.sol.RData")
constMean.truth <- list(list(rowIdx = 1:30, colIdx = 1:15))

mean(computePerf(constMean.truth, simSol$sol, amrs.hp))
median(computePerf(constMean.truth, simSol$sol, amrs.hp))
sd(computePerf(constMean.truth, simSol$sol, amrs.hp))

(computePerf(constMean.truth, simSol$sol, amrs))

corMvn.sim.test.hclust.med <- function(nSamples)
{
    set.seed(42)
    seeds <- sample.int(10000, nSamples)

    blocks <- vector(nSamples, mode = "list")
    solz <- vector(nSamples, mode = "list")

    for (it in 1:nSamples )
    {
        cat("***********************************\n")
        cat("Running simulation: ", it, "\n")
        cat("***********************************\n")
        set.seed(seeds[it])
        #curBicluster <- genBlock.iid.N(30, 15, 4, 1)
        curBicluster <- genRandomBlock.mvn(30, 15)
        curSim <- generateNormal(nrow = 300, ncol = 40, clusterOptions = 
                               list(list(x.start = 1, x.end = 30, y.start = 1, y.end = 15,
                                         data = curBicluster)))
        ps.curSim <- permuteSim(curSim)
        blocks[[it]] <- ps.curSim
        curSol <- biclusteringPar(ps.curSim$permutedMat, 8, 15)
        pss <- postSubSample(curSol)
        clusters <- post.hclust(pss)
        curSol <- list(allBcSol = list(curSol), clusters = clusters)
        curSol.reorder <- reorderBCSol(curSol, ps.curSim)

        solz[[it]] <- curSol.reorder
    }


    list(data = blocks, sol = solz)
}

withMed <- corMvn.sim.test.hclust.med(1)

noMed <- corMvn.sim.test.hclust(1)

corMvn.sim.test.hclust <- function(nSamples)
{
    set.seed(42)
    seeds <- sample.int(10000, nSamples)

    blocks <- vector(nSamples, mode = "list")
    solz <- vector(nSamples, mode = "list")

    for (it in 1:nSamples )
    {
        cat("***********************************\n")
        cat("Running simulation: ", it, "\n")
        cat("***********************************\n")
        set.seed(seeds[it])
        #curBicluster <- genBlock.iid.N(30, 15, 4, 1)
        curBicluster <- genRandomBlock.mvn(30, 15)
        curSim <- generateNormal(nrow = 300, ncol = 40, clusterOptions = 
                               list(list(x.start = 1, x.end = 30, y.start = 1, y.end = 15,
                                         data = curBicluster)))
        ps.curSim <- permuteSim(curSim)
        blocks[[it]] <- ps.curSim
        curSol <- bcMultipleClusters(ps.curSim$permutedMat, 15, 1, 8)
        curSol.reorder <- reorderBCSol(curSol, ps.curSim)
        solz[[it]] <- curSol.reorder
    }


    list(data = blocks, sol = solz)
}

twoBC <- function(nSamples)
{
    set.seed(42)
    seeds <- sample.int(10000, nSamples)

    blocks <- vector(nSamples, mode = "list")
    solz <- vector(nSamples, mode = "list")

    for (it in 1:nSamples )
    {
        cat("***********************************\n")
        cat("Running simulation: ", it, "\n")
        cat("***********************************\n")
        set.seed(seeds[it])
        addCluster <- genAdditive(50, 15)
        progCluster <- generateRandomBlock.progRows(40, 15)
        curSim <- generateNormal(nrow = 500, ncol = 40, clusterOptions = 
                               list(list(x.start = 1, x.end = 50, y.start = 1, y.end = 15,
                                         data = addCluster),
                                    list(x.start = 200, x.end = 239, y.start = 26, y.end = 40,
                                         data = progCluster))
                                 )
        ps.curSim <- permuteSim(curSim)
        blocks[[it]] <- ps.curSim
        curSol <- bcMultipleClusters(ps.curSim$permutedMat, 15, 2, 50)
        curSol.reorder <- reorderBCSol(curSol, ps.curSim)
        solz[[it]] <- curSol.reorder
    }


    list(data = blocks, sol = solz)
}
undebug(generateNormal)
time.twoBC <- system.time(twoBC.sol <- twoBC(1))
time.twoBC <- system.time(twoBC.sol <- twoBC(50))
save(twoBC.sol, file = "twoBC.sol.RData")

load("twoBC.sol.RData")
plotParSolution(twoBC.sol$sol[[1]]$allBcSol[[1]])

postSubSample.percent(twoBC.sol$sol[[1]]$allBcSol[[1]], 0.95, 0.5)


truth.twoBC <- list(list(rowIdx = 1:50, colIdx = 1:15), 
                    list(rowIdx = 200:239, colIdx = 26:40))
computePerf(truth.twoBC, twoBC.sol$sol, amrs.hp)

corMvn.sim <- function(nSamples)
{
    set.seed(42)
    seeds <- sample.int(10000, nSamples)

    blocks <- vector(nSamples, mode = "list")
    solz <- vector(nSamples, mode = "list")

    for (it in 1:nSamples )
    {
        cat("***********************************\n")
        cat("Running simulation: ", it, "\n")
        cat("***********************************\n")
        set.seed(seeds[it])
        #curBicluster <- genBlock.iid.N(30, 15, 4, 1)
        # curBicluster <- genRandomBlock.mvn(300, 30)
        curBicluster <- genRandomBlock.mvn(300, 30, .72, .94)
        curSim <- generateNormal(nrow = 1500, ncol = 100, clusterOptions = 
                               list(list(x.start = 1, x.end = 300, y.start = 1, y.end = 30,
                                         data = curBicluster)))
        ps.curSim <- permuteSim(curSim)
        blocks[[it]] <- ps.curSim
        curSol <- bcMultipleClusters(ps.curSim$permutedMat, 30, 1, 100)
        curSol.reorder <- reorderBCSol(curSol, ps.curSim)
        solz[[it]] <- curSol.reorder
    }


    list(data = blocks, sol = solz)
}


debug(genRandomBlock.mvn)
undebug(genRandomBlock.mvn)
debug(genRanPosDefMat)
mvnData <- genRandomBlock.mvn(300, 30)

mvnBlock <- generateNormal(nrow = 1500, ncol = 100, clusterOptions = 
                           list(list(x.start = 1, x.end = 300, y.start = 1, y.end = 30,
                                     data = genRandomBlock.mvn(300, 30))))
system.time(simMvn.sol.1500 <- bcMultipleClusters(mvnBlock, 30, 1, 100))
save(simMvn.sol, file = "simMvn.sol.RData")
truth.mvn <- list(list(rowIdx = 1:300, colIdx = 1:30))
amrs.hp(truth.mvn, simMvn.sol$clusters)

load("simMvn.sol.RData")
plotParSolution(simMvn.sol$allBcSol[[1]])

system.time(simMvn.sol <- corMvn.sim(1))

system.time(sim50Mvn.sol <- corMvn.sim(50))
# 231905.34
save(sim50Mvn.sol, file = "sim50Mvn.sol.RData")

load("sim50Mvn.sol.RData")
cor50Mvn.truth <- list(list(rowIdx = 1:300, colIdx = 1:30))
simMvn50.perf <- computePerf(cor50Mvn.truth, sim50Mvn.sol$sol, amrs.hp)

which(simMvn50.perf == sort(simMvn50.perf)[26])
# get a simulation near the median
nearMedMvn <- sim50Mvn.sol$sol[[41]]$allBcSol[[1]]
nearMedMvn <- lapply(nearMedMvn, function(x){
                     x$ab <- x$ab[c(1:300, sample(300:1500, 200))]
                     x
                           })
plotParSolution(nearMedMvn, "mvn5to8")
ggPlotParSolution(nearMedMvn, "mvn5to8")


simMvn.sol <- corMvn.sim(50)
save(simMvn.sol, file = "simMvn.sol.RData")


# higher cor MVN
system.time(sim50Mvn.sol.highCor <- corMvn.sim(50))

save(sim50Mvn.sol.highCor, file = "sim50Mvn.sol.highCor.RData")

load("sim50Mvn.sol.highCor.RData")

cor50Mvn.truth <- list(list(rowIdx = 1:300, colIdx = 1:30))
boxplot(computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp), ylim = c(0,1))

computePerf(cor50Mvn.truth, sim50Mvn.sol.highCor$sol, amrs.hp)

load("~/hp2.RData")
