context("Maximizing condition coefficients (Lasso)")

test_that("Regularizing correctly",
    {
        d <- .lasso_max_d(c(2, 1, -2, 30, -1, -3, 3, 3), 2.5)
        expect_equals(sum(d), 2.5)
    })

# .lasso_max_d(c(-1, -2, -10, -2, -2), 2.5)
# .lasso_max_d(c(-1, -2, -10, 20, -2, -2), 2.5)
# .lasso_max_d(c(-1, -2, -10, 20, -2, -2), 0.7)
# .lasso_max_d(c(1, 30, 30, 30, 30), 0.7)
# .lasso_max_d(c(NA, 30, 30, 30, 30), 0.7)

