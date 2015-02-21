#' Filter out rows and columns with low expression
#'
#' Filter out rows and columns with low expression
#' @param data the matrix to filter
#' @param min_exp the minimum expression in each row
#' @param max_cols the maximum number of cols that can have \code{min_exp}
#' @export
filter_low_expression <- function(data, min_exp, max_cols) {
  stopifnot(is(data, "matrix"))

  rows_kept <- apply(data, 1, function(row)
    {
      sum(row < min_exp) <= max_cols
    })

  list(data = data[rows_kept,], rows_kept = rows_kept)
}
