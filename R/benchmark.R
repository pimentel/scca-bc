const_mean_sim <- function(
    nrows_clust,
    ncols_clust,
    nrows,
    ncols,
    noise_mean = 0,
    noise_sd = 1)
{
    #generateNormal

}

#' Insert a submatrix into a matrix
#'
#' Take a matrix that is smaller and insert it into a particular position in
#' a larger matrix
#'
#' @param sub_matrix the matrix to insert
#' @param row_start the row to insert the submatrix into the background matrix
#' @param col_start the column to insert the submatrix into the background matrix
#' @param bg_matrix the background matrix

#' @return the matrix \code{bg_matrix} with \code{sub_matrix} inserted at
#' position \code{row_start} and \code{col_start}.
#' @export
insert_matrix <- function(sub_matrix, row_start, col_start, bg_matrix)
{
    if ( !is(sub_matrix, "matrix") || !is(bg_matrix, "matrix"))
        stop("'sub_matrix' and 'bg_matrix' must be of type 'matrix'")

    # if ( !is(row_start, "integer") || !is(col_start, "integer"))
    #     stop("'row_start' and 'col_start' must be of type 'integer'")

    row_end <- row_start + nrow(sub_matrix) - 1
    col_end <- col_start + ncol(sub_matrix) - 1

    if (row_end > nrow(bg_matrix))
        stop("'sub_matrix' has too many rows to fit into 'bg_matrix'")

    if (col_end > ncol(bg_matrix))
        stop("'sub_matrix' has too many cols to fit into 'bg_matrix'")

    bg_matrix[row_start:row_end, col_start:col_end] <- sub_matrix

    bg_matrix
}

#' Generate a matrix that is iid normal
#'
#' Generate a matrix that is iid normal
#'
#' @param nrows the number of rows
#' @param ncols the number of cols
#' @param mean the mean of each sample
#' @param sd the standard devian of each sample
#' @return a \code{matrix} of size nrows * ncols with iid samples from N(mean,
#' sd)
#' @export
iid_gaussian_block <- function(nrows, ncols, mean, sd)
{
    matrix(rnorm(nrows * ncols, mean = mean, sd = sd),
        nrow = nrows, ncol = ncols)
}


#' Generate a multivariate normal block
#'
#' Generate a multivariate normal block
#'
#' @param nrows the number of rows
#' @param ncols the number of cols
#' @param min the min number to pass to \link{\code{random_pos_def_mat}}
#' @param max the max number to pass to \link{\code{random_pos_def_mat}}
#' @return a multivariate normal with dimensions nrows * ncols
#' @export
mvn_block <- function(nrows = 30, ncols = 20, min = 0.5, max = 0.8)
{
    covMat <- random_pos_def_mat(nrows, min, max)
    ranSamples <- MASS::mvrnorm(ncols, mu = rep(0, nrows), Sigma = covMat)
    t(ranSamples)
}

#' Generate a random posisitive definite matrix
#'
#' Generate a random posisitive definite matrix
#'
#' @param nvars the number of variables
#' @param min the minimum number to generate the covariance from
#' @param max the maximum number to generate the covariance from
#' @return a random positive definite matrix with diagonal 1
#' @export
random_pos_def_mat <- function(nvars = 10, min = 0.5, max = 0.8)
{
    correlations <- round(runif(round((nvars^2-nvars)/2),
                                min = min, max = max), digits = 2)
    # correlations <- sample(c(1, -1), length(correlations), replace = TRUE) * correlations
    sig <- matrix(nrow = nvars, ncol = nvars)
    diag(sig) <- 1 # common variance of 1
    sig[upper.tri(sig)] <- correlations
    sig[lower.tri(sig)] <- t(sig)[lower.tri(sig)]

    # using Matrix package
    sigNear <- Matrix::nearPD(sig)

    cov2cor(as.matrix(sigNear$mat))
}


#' Permute a matrix
#'
#' Randomly permute the rows and columns of a matrix.
#'
#' @param mat any object (usually a \code{matrix}) that accepts functions \code{nrow} and \code{ncol}
#' @return a list with named elements:
#' \itemize{
#'  \item{"row_order"}{the order of the rows}
#'  \item{col_order}{the order of the columns}
#'  \item{permuted}{the permuted matrix}
#'  \item{original}{the original matrix}
#' }
#' @seealso \code{\link{translate_idx}}
permute_mat <- function(mat)
{
    row_order <- sample(nrow(mat))
    col_order <- sample(ncol(mat))

    list(row_order = row_order, col_order = col_order,
            permuted = mat[row_order, col_order],
            original = mat)
}
