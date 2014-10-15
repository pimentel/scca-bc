context("Maximizing condition coefficients (Lasso)")

test_that("Regularization",
    {
        d <- .lasso_max_d(c(2, 1, -2, 30, -1, -3, 3, 3), 2.5)
        expect_equal(sum(d), 2.5)
        expect_equal(d, c(0.0, 0.0, 0.0, 1, 0.0, 0.0, 1, 0.5))

        d <- .lasso_max_d(c(-Inf, 0, -2, -3, -0.1), 2.5)
        expect_equal(sum(d), 0.0)

        d <- .lasso_max_d(c(30, 30, 30, 30), 2.5)
        expect_equal(sum(d), 2.5)
    })

# .lasso_max_d(c(-1, -2, -10, -2, -2), 2.5)
# .lasso_max_d(c(-1, -2, -10, 20, -2, -2), 2.5)
# .lasso_max_d(c(-1, -2, -10, 20, -2, -2), 0.7)
# .lasso_max_d(c(1, 30, 30, 30, 30), 0.7)
# .lasso_max_d(c(NA, 30, 30, 30, 30), 0.7)

