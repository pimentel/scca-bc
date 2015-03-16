
context("Matrix manipulations")

test_that("Insert matrix",
    {
        bg <- matrix(0, nrow = 10, ncol = 4)
        ins <- matrix(1:8, ncol = 2, byrow = T)

        res <- insert_matrix(ins, 3, 2, bg)
        expect_equal(res[3:6, 2:3], ins)

        expect_error(insert_matrix(ins, 3, 4, bg))

        expect_error(insert_matrix(ins, 10, 2, bg))
    })

test_that("permuting a matrix",
  {
    mat <- matrix(1:100, nrow = 20)
    set.seed(42)
    pm <- permute_mat(mat)
    expect_false(identical(mat, pm$permuted))
    expect_equal(unpermute_mat(pm$permuted, pm), mat)
  })
