library("devtools")
library("sccab")

gen_high_cor <- function(n_samps) {
  sims <- lapply(1:n_samps, function(i)
    {
      bg_matrix <- iid_gaussian_block(nrows = 1500, ncols = 100, mean = 0, sd = 1)
      cur_mvn <- mvn_block(nrows = 300, ncols = 30, min = 0.72, max = 0.94)
      sim_data <- insert_matrix(cur_mvn, 1, 1, bg_matrix)

      sim_data
    })

  structure(list(sims = sims, truth = list(list(rowIdx = 1:300, colIdx = 1:30))),
    class = c("sccab_high_cor", "sccab_sim"))
}

p_res <- pps(res, "hclust", "hclust")
jaccard_idx_matrix(p_res, truth)

test("~/dev/sccab")

document("~/dev/sccab")
install("~/dev/sccab")

x <- gen_high_cor(5)

params <- sccab_params(d_upr = 30, d_lwr = 0, n_samp = 10)
s_res <- sccab_sims(x, "junk", params)

p_res <- lapply(s_res, function(s) pps(s, "hclust", "hclust")) %>%
  lapply(function(s) jaccard_idx_matrix(s, x$truth))
