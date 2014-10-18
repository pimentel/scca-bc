vScore <- function(set1, set2)
{
    length(intersect(set1, set2)) / length(union(set1, set2))
}

vScore2 <- function(...)
{
    length(intersect(...)) / length(union(...))
}

# TODO: Replace the regular distance in maximizeOneSplit with pDiff and see how
# affects convergence
p_diff <- function(old, new)
{
    num <- dist(rbind(old, new))[1]
    denom <- sqrt(sum(old^2))
    num / denom
}

mean_relative_tol <- function(old, new)
{
    mean(abs(1 - new / old ))
}

mean_absolute_tol <- function(old, new)
{
    mean(abs(new - old))
}

#' @export
jaccard_idx_matrix <- function(truthList, predList)
{
    expandCells <- function(aList)
    {
        aCat <- expand.grid(aList$rowIdx, aList$colIdx)
        paste(aCat[,1], aCat[,2], sep = ".")
    }
    mean(sapply(truthList, function(truth)
           {
               truthCat <- expandCells(truth)
               max(sapply(predList, function(prediction)
                      {
                          predCat <- expandCells(prediction)
                          numerator <- length(intersect(truthCat, predCat)) 
                          denominator <- length(union(truthCat, predCat))
                          numerator / denominator
                      }))
           }))
}

