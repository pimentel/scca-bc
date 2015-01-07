#' Plot the parallel solution
#'
#' Plot a solution from either bcSubSampleSerial or bcSubSamplePar
#'
#' @param parSol the solution
#' @param clust if you've post-process clustered the clustered object denoting
#' which rows/columns exist in the bicluster
#' @param fBase a string for the file name
#' @param rowNames if not NA, the rownames to use for expression features
#' @param colNames if not NA, a character vector of column (condition) names
#' @return A list with two ggplot objects
#' @export
ggPlotSSSolution <- function(parSol, clust = NA, fBase = NA, rowNames = NA,
                              colNames = NA, re_order = NULL)
{
    par(ask = FALSE)

    A <- parSol$AB
    D <- parSol$D

    if (!is.na(clust)[1])
    {
        # TODO: make sure valid cluster (i.e. has rowIdx, colIdx)
        rowOrder <- c(clust$rowIdx, setdiff(1:nrow(A), clust$rowIdx))
        colOrder <- c(clust$colIdx, setdiff(1:nrow(D), clust$colIdx))
        A <- A[rowOrder,]
        D <- D[colOrder,]
    }
    dRot <- 0
    dxSize <- 14
    fSize <- 1.5
    if (!is.na(rowNames)[1])
    {
        rownames(A) <- rowNames
    }
    if (!is.na(colNames)[1])
    {
        rownames(D) <- colNames
        dRot <- 90
        dxSize <- 10
        fSize = 1
    }
    meltA <- melt(t(A))
    colnames(meltA) <- c("x", "y", "value")
    breaksA <- round(seq(0, max(meltA$value, na.rm = T), length.out = 10), 3)
    p <- ggplot(meltA, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradient(low = "firebrick2", high = "yellow", 
                                               guide = guide_legend(title = "Coefficient", reverse = T), 
                                               breaks = breaksA) 

    p <- p + xlab("Iteration") + ylab("Feature") + theme_bw() + 
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=14), 
              axis.title=element_text(size=15), 
              legend.text=element_text(size=14), 
              legend.title=element_text(size=14))
    res <- list(featurePlot = p)

    print(p)
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        ggsave(paste("../img/gg", fBase, "Feature.pdf", sep = ""), 
               width = 6.62, height = 12.8)
    meltD <- melt(D)
    colnames(meltD) <- c("x", "y", "value")
    breaksD <- round(seq(0, max(meltD$value, na.rm = T), length.out = 10), 1)
    p <- ggplot(meltD, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradient(low = "firebrick2", high = "yellow", 
                                               guide = guide_legend(title = "Coefficient", reverse = T), 
                                               breaks = breaksD) 
    p <- p + xlab("Condition") + ylab("Iteration") + theme_bw() + theme(axis.text.x=element_text(size=dxSize,
                                                                                                 angle = dRot),
                                                    axis.text.y=element_text(size=14), 
                                                    axis.title =element_text(size=15),
                                                    legend.text=element_text(size=14),
                                                    legend.title=element_text(size=14))
    print(p)
    readline("Next [Enter]\t")
    if (!is.na(fBase))
        ggsave(paste("../img/gg", fBase, "Cond.pdf", sep = ""),
               width = 6.62, height = 12.8)

    res$condPlot <- p
    # Saving 6.62 x 12.8 in image
    invisible(res)
}




condSizePlot <- function(pssSols, dat, range = NULL, truth = NULL)
{
    if (is.null(range) || length(pssSols) != length(range))
        stop("Range must be a non-null value of the lambda values used, 
             the same length of pssSols")

    sumDf <- ldply(pssSols, function(sol)
                   {
                       corDist <- cor(t(dat[sol$rowIdx, sol$colIdx]))
                       corDist <- abs(corDist[upper.tri(corDist)])

                       curAmrs <- NA

                       if (!is.null(truth))
                           curAmrs <- amrs.hp(truth, list(sol)) 


                       return(c(mean(corDist), median(corDist), sd(corDist),
                                    length(sol$rowIdx), length(sol$colIdx), 
                                    curAmrs))
                       # return (c(mean(corDist), median(corDist), sd(corDist), 
                       #           length(sol$rowIdx), length(sol$colIdx)
                       #           ))
                   })

    pairAmrs <- sapply(adjacentPairs(1, length(pssSols)), 
                       function(x) amrs.hp(list(pssSols[[x[1]]]), list(pssSols[[x[2]]])))
    pairAmrs <- c(NA, pairAmrs)
    colV <- sapply(adjacentPairs(1, length(pssSols)), 
                   function(x) vScore(pssSols[[x[1]]]$colIdx, pssSols[[x[2]]]$colIdx))
    colV <- c(NA, colV)
    rowV <- sapply(adjacentPairs(1, length(pssSols)), 
                   function(x) vScore(pssSols[[x[1]]]$rowIdx, pssSols[[x[2]]]$rowIdx))
    rowV <- c(NA, rowV)
    vIdx <- c(NA, (range + 0.5)[1:(length(range))-1])


    rownames(sumDf) <- range
    sumDf <- cbind(range, sumDf, vIdx, pairAmrs, colV, rowV)
    colnames(sumDf) <- c("lambda", "mean", "median", "sd", "nRows", "nCols", "amrs",
                         "vIdx", "pairAmrs", "colV", "rowV")

    meanPlot <- ggplot(sumDf, aes(x = lambda)) + 
        geom_line(aes(y = mean, colour = "mean")) + 
        geom_point(aes(y = mean, colour = "mean")) + 
        geom_line(aes(y = amrs, colour = "amrs"), na.rm = T) + 
        geom_point(aes(y = amrs, colour = "amrs"), na.rm = T) + 
        geom_line(aes(y = sd, colour = "sd")) + 
        geom_point(aes(y = sd, colour = "sd")) + 
        # geom_line(aes(y = median, colour = "median")) + 
        # geom_point(aes(y = median, colour = "median")) + 
        geom_line(aes(x = vIdx, y = pairAmrs, colour = "pairAmrs"), na.rm = T) + 
        geom_point(aes(x = vIdx, y = pairAmrs, colour = "pairAmrs"), na.rm = T) +
        # geom_line(aes(x = vIdx, y = colV, colour = "colV"), na.rm = T) + 
        # geom_point(aes(x = vIdx, y = colV, colour = "colV"), na.rm = T) +
        # geom_line(aes(x = vIdx, y = rowV, colour = "rowV"), na.rm = T) + 
        # geom_point(aes(x = vIdx, y = rowV, colour = "rowV"), na.rm = T) +
        scale_y_continuous(limits = c(0, 1))

    
    rowPlot <- ggplot(sumDf, aes(lambda, nRows)) + 
        geom_line(aes(colour = "nRows")) + geom_point(aes(colour = "nRows"))
    colPlot <- ggplot(sumDf, aes(x = lambda, y = nCols)) + 
        geom_line(aes(colour = "nCols")) + geom_point(aes(colour = "nCols")) +
        geom_line(aes(y = lambda, colour = "lambda"))


    amrsPlot <- NA
    if (!is.null(truth))
        amrsPlot <- ggplot(sumDf, aes(lambda, amrs)) + geom_line() + 
            geom_point() + 
            scale_y_continuous(limits = c(0,1))

    return(list(summaryDf = sumDf, 
                meanPlot = meanPlot,
                rowPlot = rowPlot,
                colPlot = colPlot,
                amrsPlot = amrsPlot))
}




plotHeatmap <- function(data, clust, rows = T, cols = F)
{
    if (rows)
        rows <- NULL
    else
        rows <- NA
    if (cols)
        cols <- NULL
    else
        cols <- NA

    # heatmap(data[clust$rowIdx, clust$colIdx], Rowv = rows, Colv = cols, 
    #         col = redgreen(256), keep.dendro = F, scale = "none")

    # heatmap.2(data[clust$rowIdx, clust$colIdx], Rowv = rows, Colv = cols, 
    #         col = redgreen(256), dendrogram = "both", scale = "none")
    heatmap.2(data[clust$rowIdx, clust$colIdx], col = redgreen(256), trace = "none")
}

plotHeatmap2 <- function(data, rows = T, cols = F)
{
    if (rows)
        rows <- NULL
    else
        rows <- NA
    if (cols)
        cols <- NULL
    else
        cols <- NA

    heatmap(data, Rowv = rows, Colv = cols, col = redgreen(256))
}


corFigure <- function(sol, data, truth = FALSE)
{
    corClust <- cor(t(data[sol$rowIdx, sol$colIdx]))
    corAll <- cbind(corClust[upper.tri(corClust)], "Bicluster")

    ranRow <- sample.int(n = nrow(data), size = length(sol$rowIdx))
    ranRowCor <- cor(t(data[ranRow, sol$colIdx]))
    corAll <- rbind(corAll, cbind(ranRowCor[upper.tri(ranRowCor)], "Random rows"))


    ranRow <- sample.int(n = nrow(data), size = length(sol$rowIdx))
    ranCol <- sample.int(n = ncol(data), size = length(sol$colIdx))
    bothRanCor <- cor(t(data[ranRow, ranCol]))
    corAll <- rbind(corAll, cbind(bothRanCor[upper.tri(bothRanCor)], "Both random"))

    allCond <- cor(t(data[sol$rowIdx,]))
    corAll <- rbind(corAll, cbind(allCond[upper.tri(allCond)], "All conditions"))

    if (truth)
    {
        corTruth <- cor(t(data[1:300, 1:30]))
        corAll <- rbind(corAll, cbind(corTruth[upper.tri(corTruth)], "Truth"))
    }

    corAll <- as.data.frame(corAll, stringsAsFactors = F)
    corAll[,1] <- as.numeric(corAll[,1])
    colnames(corAll) <- c("Correlation", "Conditions")
    
    # ggplot(corAll, aes(x = Conditions, y = abs(Correlation), colour = Conditions)) + 
    # ggplot(corAll, aes(x = Conditions, y = abs(Correlation), colour = Conditions)) + 
    ggplot(corAll, aes(x = abs(Correlation), y = ..density.., colour = Conditions)) + 
        # geom_density(aes(fill = Conditions), alpha = 0.5) + 
        geom_histogram(aes(fill = Conditions), position = "dodge") + 
        # geom_boxplot(aes(fill = Conditions), alpha = 0.5) + 
        scale_fill_manual(values=cbbPalette) + 
        scale_colour_manual(values=cbbPalette)
}

bcGeneNames <- function(data, sol = NULL)
{
    if (is.null(sol))
        return(rownames(data))
    rownames(data[sol$rowIdx,])
}

genePlot <- function(data, sol, genes, ordering = NULL)
{
    # sol$colIdx <- colnames(data)[sol$colIdx]
    bg <- cbind(x = 1:ncol(data), melt(t(data[genes,])))
    colnames(bg) <- c('x', 'condition', 'gene', 'value')
    fg <- cbind(x = sol$colIdx, melt(t(data[genes, sol$colIdx])))
    colnames(fg) <- c('x', 'condition', 'gene', 'value')
    plt <- ggplot(bg, aes(x = x, y = value, color = gene)) + geom_line()
    plt <- plt + geom_point(aes(x = x, y = value), data = fg,
                            size = 3, shape = 21, fill = 'black') 
    plt
}

#' Plot the expression only in the bicluster
#'
#' Given some solution from a post-process method, plot only those rows and
#' columns
#'
#' @param exMat the expression matrix
#' @param clust the post-process clustering output
#' @return a gpplot object
#' @export
plotClusterExpression <- function(exMat, clust)
{
    ggPlotExpression(dat[clust$rowIdx, clust$colIdx])
}

#' Heatmap of expression
#'
#' Plot all of the points in an expression matrix
#'
#' @param exMat the expression matrix
#' @param cluster either the cluster results to group the rows/columns in the
#' biclustering, or NULL if there are no results

#' @param clustRows if TRUE, cluster the rows by hierarchical clustering.
#' @param clustCols if TRUE, cluster the columns by hierarchical clustering.
#' @param rowNames if TRUE, print the row names on the plot
#' @param colNames if TRUE, print the column names on the plot
#' @return a ggplot object
#' @export
ggPlotExpression <- function(exMat, cluster = NULL, clustRows = T, clustCols = T,
                             rowNames = F, colNames = T)
{
    if (class(exMat) != 'matrix')
    {
        exMat <- as.matrix(exMat)
        stopifnot(class(exMat) == 'matrix')
    }
    if (!is.null(cluster))
        exMat <- exMat[cluster$rowIdx, cluster$colIdx]
    rowOrder <- 1:nrow(exMat)
    colOrder <- 1:ncol(exMat)
    if (clustRows)
        rowOrder <- orderByDendrogram(exMat)
    if (clustCols)
        colOrder <- orderByDendrogram(t(exMat))
    exMat <- exMat[rowOrder, colOrder]
    meltMat <- reshape2::melt(exMat, varnames = c("x", "y"))
    breaksM <- round(seq(min(meltMat$value, na.rm = T), max(meltMat$value, na.rm = T), 
                         length.out = 10), 3)
    print(rownames(exMat))
    if (is.null(colnames(exMat)))
        colnames(exMat) <- 1:ncol(exMat)
    meltMat$y <- factor(meltMat$y, levels = colnames(exMat))
    p <- ggplot(meltMat, aes(x, y, fill = value))
    p <- p + geom_tile() + scale_fill_gradientn(colours = redgreen(20),
                                                guide = guide_legend(title = "Expression", 
                                                                     reverse = T, size = 14)) 
    p <- p + xlab("Gene") + ylab("Condition") + theme_bw() + theme(legend.text = element_text(size = 14),
                                                                   legend.title = element_text(size = 14),
                                                                   axis.title=element_text(size=15))
    if (rowNames)
        p <- p + theme(axis.text.x=element_text(angle = 90, size=14))
    else
        p <- p + theme(axis.text.x=element_text(size=0))

    if (colNames)
        p <- p + theme(axis.text.y=element_text(size=14))
    else
        p <- p + theme(axis.text.y=element_text(size=0))


    list(plot = p, rowOrder = rowOrder, colOrder = colOrder)
}

#' Order by dendrogram
#'
#' @param mat a matrix where the rows are observations and the columns are different dimensions on the matrix
#' @return a vector of label orderings
#' @export
orderByDendrogram <- function(mat)
{
    hc <- hclust(dist(mat))
    dc <- as.dendrogram(hc)
    order.dendrogram(dc)
}
