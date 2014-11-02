#' Optimizing the biclustering function
#'
#' One iteration of biclustering maximization. Takes the current values of
#' d, and if successful returns one iteration.
#'
#' @param X matrix with rows representing genes and columns representing conditions (n x k)
#' @param Y matrix with rows representing genes and columns representing conditions (n x k)
#' @param d_init vector of dimension k
#' @param ab_lam regularization on expression features
#' @param ab_lam regularization on expression features
#' @return A list containing maximum values of a, b, d
#' @export
biClustMax.optim <- function(X, Y, d_init, ab_lam, d_lwr, d_upr, d_max_fun)
{
    # First maximize a and b using SCCA
    res <- features_max_fscca(X, Y, d_init, ab_lam, ab_lam)

    d_optim <- d_max_fun(X, Y, res$a, res$b, d_upr)

    # print(res$lambda)

    list(a = res$a, b = res$b, d = d_optim, q, sccaLam = res$lambda)
}

#' Maximize one split of the SCCAB algorithm
maximizeOneSplit <- function(exp_mat, params, epsA = 0.001, epsB = 0.001,
    epsD = 0.01, maxIt = 100)
{
    splitIdx <- split_evenly(nrow(exp_mat))
    X <- exp_mat[splitIdx[[1]], ]
    Y <- exp_mat[splitIdx[[2]], ]

    # randomly initialize d
    # FIXME: If there are lots of conditions, possible that won't be able to
    # find d that satisfy the constraint (if lam is sufficiently small. Come up
    # with different way to randomly assign values
    # d <- 0
    # minUnif <- 0.05
    # maxUnif <- 1

    # the first iteration, simply compute SCCA and get the highest ranked
    # conditions

    d <- rep.int(1, ncol(X))

    # repeat {
    #     d <- runif(ncol(X), max = maxUnif)
    #     if (sum(d) <= lam)
    #         break
    #     else
    #         maxUnif <- min(minUnif + 0.01, maxUnif - 0.05)
    # }

    # NB: values of a and b are not important currently since evaluated after
    # d is set.  If ever change the order, fix this.
    a <- b <- runif(nrow(X), min = -1, max = 1)
    it <- 1
    curSol <- list()
    cat(sprintf("% 10s% 7s% 7s\n", "a", "b", "d"))
    repeat {
        # cat("\tOne split iteration: ", it, "\n")
        curSol <- biClustMax.optim(X, Y, d, params$ab_lam, params$d_lwr,
            params$d_upr, params$optim_fun)

        a_dist <- dist(rbind(curSol$a, a))[1]
        b_dist <- dist(rbind(curSol$b, b))[1]
        # d_dist <- dist(rbind(curSol$d, d))[1]
        d_dist <- mean_absolute_tol(d, curSol$d)

        cur_val <- (t(a) %*% (X %*% diag(d))) %*% t(t(b) %*% (Y %*% diag(d)))

        cat(sprintf("\t%.4f\t%.4f\t%.4f\t%.4f\n",
                a_dist, b_dist, d_dist, cur_val))

        if (a_dist < epsA && b_dist < epsB && d_dist < epsD && it >= 3)
        {
            cat("\tConverged. Number of iterations: ", it, "\n")
            break
        }

        a <- curSol$a
        b <- curSol$b
        d <- curSol$d

        it <- it + 1
        if (it > maxIt)
        {
            cat("Warning: exceeded maximum number of iterations (", maxIt, ")\n")
            break
        }
    }

    # NB: If there is something funky with the distribution of a and b,
    # check here... though this should work correctly
    ab <- rep.int(0, length(b)*2)
    ab[splitIdx[[1]]] <- curSol$a
    ab[splitIdx[[2]]] <- curSol$b

    list(a = a, b = b, ab = ab, d = round(curSol$d, digits = 4),
        splitIdx = splitIdx, sccaLam = curSol$sccaLam)
}

#' SCCAB(iclustering) (with all data)
#'
#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(cores = N).
#'
#' @param exp_mat matrix with features on rows and conditions on columns
#' @param params a sccab_params object from sccab_params()
#' @return a sccab_res object
#' @export
sccab <- function(exp_mat, params)
{
    if (!is.matrix(exp_mat))
        stop("sccab requires a matrix")

    if (class(params) != "sccab_params")
        stop("params must be a 'sccab_params' object. Please run sccab_params()")

    p <- params

    check_lams(p$d_lwr, p$d_upr, ncol(exp_mat))

    p$apply_fun(1:p$n_samp,
        function(it) {
            cat("Biclustering iteration: ", it, "\n")
            maximizeOneSplit(exp_mat, p)
        })
}

#' SCCAB(iclustering) with subsampling
#'
#' Given a set of features, will find the most "dominant" bicluster. Works in
#' parallel using library "multicore." To set the number of cores used, set
#' options(mc.cores = N).
#'
#' @param exp_mat matrix with features on rows and conditions on columns
#' @param params a sccab_params object from sccab_params()
#' @return a sccab_res object
#' @export
sccab_subsample <- function(exp_mat, params)
{
    if (!is.matrix(exp_mat))
        stop("sccab requires a matrix")

    p <- params

    if (is.na(p$prop))
        stop("Must set a sampling proportion (prop) in sccab_params()")

    check_lams(p$d_lwr, p$d_upr, ncol(exp_mat))

    nRowsSample <- round(p$prop * nrow(exp_mat))
    if (nRowsSample %% 2)
        nRowsSample <- nRowsSample + 1

    cat("Sampling ", nRowsSample, " features\n")
    p$apply_fun(1:p$n_samp, function(it) {
        cat("Biclustering iteration: ", it, "\n")
        sampIdx <- sample.int(nrow(exp_mat), size = nRowsSample)
        abSol <- rep.int(NA, nrow(exp_mat))
        curSol <- maximizeOneSplit(exp_mat[sampIdx,], p)

        abSol[sampIdx] <- curSol$ab
        d <- curSol$d

        list(ab = abSol, d = d, sccaLam = curSol$sccaLam)
    })
}
