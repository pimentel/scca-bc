check_lams <- function(lam_lwr, lam_upr, n_cols)
{
    cur_values <- sprintf("Current values: (%.4f, %.4f)", lam_lwr, lam_upr)

    if (lam_lwr < 0)
        stop(sprintf("lam_lwr must be >= 0. %s", cur_values))

    if (!(lam_lwr - 1 <= lam_upr))
        stop(sprintf("lam_upr must be greater than or equal to lam_lwr - 1.
                %s", cur_values))

    if (lam_upr >= n_cols)
        warning("lam_upr is >= n_cols.")
}

#' Randomly split a set of indices
#'
#' Randomly splits indices if there are an even number. If there is an odd
#' number, it will stop().
#'
#' @param n the number of indices
#' @return a list with two disjoint sets of indices from 1 to n.
split_evenly <- function(n)
{
    if ( n %% 2 != 0)
        stop("Not an even number of samples!")

    sets <- split(sample.int(n), 1:2)

    sets
}

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
