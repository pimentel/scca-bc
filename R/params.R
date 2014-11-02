detect_d_maximizer <- function(opt)
{
    if (class(opt) != "character" && class(opt) != "function")
        stop("Variable 'opt' must be a character (\"lasso\", \"timeseries\") " +
            " or function")

    if (class(opt) == "function")
        return(opt)

    switch (
        opt,
        lasso = {
            return(lasso_max_d)
        },
        timeseries = {
            return(timeseries_max_d)
        },
        stop("Unrecognized 'opt' type. Please specify \"lasso\"," +
            " \"timeseries\", or pass a function")
        )
}

#' Set SCCAB(biclustering) parameters
#'
#' Call this function before calling any sccab function. The result is used in
#' the sccab functions to define behavior.
#'
#' @param d_upr regularization parameter for conditions (maximum number of conditions allows in a bicluster)
#' @param n_samp integer denoting the number of permutations to perform
#' @param d_lwr regularization parameter for conditions (minimum number of conditions allows in a bicluster)
#' @param ab_lam vector of regularization parameters for expression features.
#' @param optim_fun "lasso" or "timeseries" or a user-supplied function. NOTE: there are no checks for validity of the user-supplied function.
#' @param prop proportion to subsample
#' @param parallel if TRUE use parallel::mclapply, else use lapply
#' @param verbose if TRUE print messages
#' @return A sccab_params object to be used in sccab or sccab_subsample
#' @export
sccab_params <- function(d_upr, n_samp = 100, d_lwr = 3, ab_lam = c(0.5,3, 10),
    optim_fun = "lasso", prop = NA, parallel = FALSE, verbose = TRUE)
{
    res <- list()
    if (d_lwr < 0)
        stop(sprintf("d_lwr must be at least zero: (%.2f)\n", d_lwr))

    if (d_lwr > d_upr)
        stop(sprintf("d_lwr must be less than d_upr: (%.2f, %.2f)\n",
                d_lwr, d_upr))

    res$d_lwr <- d_lwr
    res$d_upr <- d_upr

    if (class(ab_lam) != "numeric" && class(ab_lam) != "integer")
        stop("ab_lam must be integers or numerics")

    if (any(ab_lam <= 0))
        stop("all ab_lam must be greater than or equal to zero")

    res$ab_lam <- ab_lam

    n_samp <- as.integer(n_samp)

    if (n_samp < 1)
        stop(sprintf("n_samp must be at least one: (%d)\n", n_samp))

    res$n_samp <- n_samp

    res$optim_fun <- detect_d_maximizer(optim_fun)

    if (!is.na(prop))
    {
        if (class(prop) != "numeric")
            stop("prop must be a numeric")

        if (prop == 1.0)
            warning("For prop == 1, use sccab")

        if (prop > 1 | prop < 0.01)
            stop("Invalid range for prop [0.01, 1)")
    }
    res$prop <- prop

    if (class(parallel) != "logical")
        stop("parallel must be either TRUE or FALSE")

    if (parallel)
        res$apply_fun <- parallel::mclapply
    else
        res$apply_fun <- lapply

    if (class(verbose) != "logical")
        stop("verbose must be either TRUE or FALSE")

    res$verbose <- verbose

    class(res) <- "sccab_params"

    res
}
