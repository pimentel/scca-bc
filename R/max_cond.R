lasso_max_d <- function(X, Y, a, b, s)
{

}

#' Maximize the d vector in biclustering
#'
#' Find the solution to the quadratic program:
#' f(d) = d^2 * q
#' s.t. |d| <= s, 0 <= d_i <= 1
#'
#' @param q a vector of length k, where k is the number of conditions
#' @param s a regularization parameter such that s > 0 and t(d) *%* q %* d  <= s
#' @return a valid maximum  
.lasso_max_d <- function(q, s)
{
    # XXX: Not checking if s > 0... a program above should check

    d <- rep.int(0, length(q))

    n_gt0 <- length(which(q > 0))

    # no terms are positive, thus the only way to maximize is to choose all
    # values to be zero
    if (n_gt0 == 0 )
        return(d)

    q_order <- order(q, decreasing = T)
    cur_sum <- 0
    i <- 1
    for (i in 1:n_gt0)
    {
        if (cur_sum + 1 <= s)
        {
            cur_sum <- cur_sum + 1
            d[q_order[i]] <- 1
        }
        else {
            d[q_order[i]] <- s - cur_sum
            break
        }
    }

    d
}
