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
