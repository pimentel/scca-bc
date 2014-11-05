#' Detect the conditions maximizer to use
#'
#' Given a string or a function, return a function to maximize the condition
#' coefficients.
#'
#' @param opt a \code{string} either "lasso" or "timeseries". If
#' a \code{function}, simply use the \code{function}. Does not check the return
#' type if \code{opt} is a \code{function}.
#' @return a function to maximize the condition coefficients
detect_d_maximizer <- function(opt)
{
    if ( !is(opt, "character") && !is(opt, "function") )
        stop("Variable 'opt' must be a character (\"lasso\", \"timeseries\") " +
            " or function")

    if (is(opt, "function"))
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

    if ( !is(ab_lam, "numeric") && !is(ab_lam, "integer"))
        stop("ab_lam must be integers or numerics")

    if (any(ab_lam <= 0))
        stop("all ab_lam must be greater than or equal to zero")

    res$ab_lam <- ab_lam

    n_samp <- as.integer(n_samp)

    if (n_samp < 1)
        stop(sprintf("n_samp must be at least one: (%d)\n", n_samp))

    res$n_samp <- n_samp

    res$optim_fun <- detect_d_maximizer(optim_fun)

    if (is(optim_fun, "function"))
        res$optim_str <- "User-supplied"
    else
        res$optim_str <- optim_fun

    if (!is.na(prop))
    {
        if ( !is(prop, "numeric") )
            stop("prop must be a numeric")

        if (prop == 1.0)
            warning("For prop == 1, use sccab")

        if (prop > 1 | prop < 0.01)
            stop("Invalid range for prop [0.01, 1)")

        if (prop < 0.55)
            warning("Recommend prop >= 0.55")
    }
    res$prop <- prop

    if ( !is(parallel, "logical") )
        stop("parallel must be either TRUE or FALSE")

    if (parallel)
        res$apply_fun <- parallel::mclapply
    else
        res$apply_fun <- lapply

    if ( !is(verbose, "logical") )
        stop("verbose must be either TRUE or FALSE")

    res$verbose <- verbose

    class(res) <- "sccab_params"

    res
}

#' Print sccab_params object
#'
#' Print \code{sccab_params} object and returns it \emph{invisibly} (via
#' \code{\link{invisible}(x)}).
#' @param \code{obj} an object of type \code{sccab_params}
#' @return the original \code{obj}
#' @export
print.sccab_params <- function(obj)
{
    cat("sccab_params\n")
    cat("----------------------------------------\n")
    cat(sprintf("%10s:\t(%.2f, %.2f)\n", "d", obj$d_lwr, obj$d_upr))
    cat(sprintf("%10s:\t(%s)\n", "ab", paste(obj$ab_lam, collapse = ", ")))
    if (!is.na(obj$prop))
        cat(sprintf("%10s:\t%.2f\n", "prop", obj$prop))
    cat(sprintf("%10s:\t%d\n", "n_samp", obj$n_samp))
    cat(sprintf("%10s:\t%s\n", "optim", obj$optim_str))

    if (identical(obj$apply_fun, lapply))
        cat(sprintf("%10s:\t%s\n", "parallel", FALSE))
    else
        cat(sprintf("%10s:\t%s\n", "n-cores", getOption("mc.cores")))

    cat(sprintf("%10s:\t%s\n", "verbose", obj$verbose))

    invisible(obj)
}
