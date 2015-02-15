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
