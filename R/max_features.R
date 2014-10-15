#' Maximize the features using fscca
#'
#' Wrapper to maximize features using fscca
#' @param X matrix of size n x p
#' @param Y matrix of size n x p
#' @param d vector of length p denoting inclusion
#' @param lam_a a vector of parameters to regularize a
#' @param lam_b a vector of parameters to regularize b
#' @return a list with 'a' and 'b' solutions from SCCA
features_max_fscca <- function(X, Y, d, lam_a, lam_b)
{
    xTilde <- t(X %*% diag(d))
    yTilde <- t(Y %*% diag(d))

    res <- fscca(xTilde, yTilde, "lasso", "lasso", lam_a, lam_b)

    return(list(a = res$A[,1], b = res$B[,1]))
}
