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

context("Maximizing feature coefficients (SCCA)")

test_that("Correct dimension",
    {
        n <- 300
        k <- 30
        X <- matrix(rnorm(n * k), ncol = k)
        Y <- matrix(rnorm(n * k), ncol = k)
        d <- runif(k)
        lam <- 1:3

        capture.output(res <- features_max_fscca(X, Y, d, lam, lam))

        expect_equal(length(res$a), n)
        expect_equal(length(res$b), n)
    })
