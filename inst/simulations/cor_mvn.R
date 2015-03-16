library("devtools")
library("sccab")

gen_high_cor <- function(n_samps) {
  sims <- lapply(1:n_samps, function(i)
    {
      bg_matrix <- iid_gaussian_block(nrows = 1500, ncols = 100, mean = 0, sd = 1)
      cur_mvn <- mvn_block(nrows = 300, ncols = 30, min = 0.72, max = 0.94)
      sim_data <- insert_matrix(cur_mvn, 1, 1, bg_matrix)

      permute_mat(sim_data)
    })

  structure(list(sims = sims, truth = list(list(rowIdx = 1:300, colIdx = 1:30))),
    class = c("sccab_high_cor", "sccab_sim"))
}

test("~/dev/sccab")

document("~/dev/sccab")
install("~/dev/sccab")

# 'hc' stands for high_cor
set.seed(42)
hc_sims <- gen_high_cor(100)

hc_params <- sccab_params(d_upr = 30, d_lwr = 0, n_samp = 50, parallel = TRUE)
hc_sccab_res <- sccab_sims(hc_sims, hc_params)

hc_pps_hclust <- mclapply(hc_sccab_res, function(s) pps(s, "hclust", "hclust"))
hc_pps_hclust <- lapply(seq_along(hc_pps_hclust),
  function(i) translate_solution(hc_pps_hclust[[i]], hc_sims$sims[[i]]))
hc_hclust_jidx <- unlist(lapply(hc_pps_hclust, function(s) jaccard_idx_matrix(s, hc_sims$truth)))


document("~/sccab")
install("~/sccab")
