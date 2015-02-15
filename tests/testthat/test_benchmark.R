
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
