# TODO: include minimum condition size to be 3
# TODO: run this for 20 - 50
# TODO: For selecting the "optimal" condition, look at average correlation plot
# (should increase to a flat region) and also look at the number of conditions
# selected. Should also flatten out

# TODO: instead of minLam and maxLam take a vector of lams
optConditionSize <- function(df, minLam, maxLam, cutoff = 0.6, 
                             bcMethod = bcSubSamplePar, 
                             postProcess = postSubSample.pca, ...)
{
    D <- as.data.frame(matrix(0, nrow = maxLam - minLam + 1, ncol = ncol(df)))
    rownames(D) <- minLam:maxLam
    sols <- list()
    it <- 1
    for (l in minLam:maxLam)
    {
        cat("Doing optimization for lambda = ", l, "\n")

        # temporarily supress output
        # sink("/dev/null")
        compTime <- system.time({
            curSol <- bcMethod(df, lam = l, lam.lwr = 3)
            sols[[it]] <- curSol
            if (!is.null(postProcess))
            {
                curClust <- postProcess(curSol, ...)
                D[it, curClust$colIdx] <- 1
            }
            it <- it + 1
        })
        # sink()
        print(compTime)
        cat(length(curClust$colIdx), "\n")
    }

    return(list(sol = sols, D = D))
}
