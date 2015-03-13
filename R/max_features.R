#' Maximize the features using fscca
#'
#' Wrapper to maximize features using fscca
#' @param X matrix of size n x p
#' @param Y matrix of size n x p
#' @param d vector of length p denoting inclusion
#' @param lam_a a vector of parameters to regularize a
#' @param lam_b a vector of parameters to regularize b
#' @return a list with 'a' and 'b' solutions from SCCA
#' @export
features_max_fscca <- function(X, Y, d, lam_a, lam_b)
{
    # TODO: can we just vectorize this?
    xTilde <- t(X %*% diag(d))
    yTilde <- t(Y %*% diag(d))

    capture.output(res <- fscca::fscca(xTilde, yTilde, "lasso", "lasso", lam_a, lam_b))
    # res <- fscca::fscca(xTilde, yTilde, "lasso", "lasso", lam_a, lam_b)

    return(list(a = res$A[,1], b = res$B[,1], lambda = res$lambda))
}

# the old functionality
    # if (verbose)
    #     cat("Using penalty: ", pen, "\n")

    # sccaRes <- scca(as.matrix(xTilde), as.matrix(yTilde), 
    #                 nc = 1, 
    #                 penalty = pen,
    #                 lamx = lamx,
    #                 lamy = lamx,
    #                 tuning = cv,
    #                 center = TRUE, scale = TRUE)
    # a <- as.numeric(sccaRes$A)
    # b <- as.numeric(sccaRes$B)

    # if (verbose)
    #     cat("Computing SCCA component\n")


