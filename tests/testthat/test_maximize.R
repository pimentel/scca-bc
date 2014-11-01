context("Max condition coef (Lasso)")

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

context("Max feature coef (SCCA)")

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

context("Timeseries max d")
test_that("Maximize",
    {
        q <- c(0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1)
        d <- .timeseries_max_d(q, 3)
        expect_equal(d, c(rep(0, 4), rep(1, 3), rep(0, 4)))
    })

context("Misc supp functions")

test_that("Detect d maximizer",
    {
        expect_equal(detect_d_maximizer("lasso"), lasso_max_d)
        expect_equal(detect_d_maximizer("timeseries"), timeseries_max_d)

        expect_equal(detect_d_maximizer(lasso_max_d), lasso_max_d)
        expect_equal(detect_d_maximizer(timeseries_max_d), timeseries_max_d)

        expect_error(detect_d_maximizer("lasso1"))
        expect_error(detect_d_maximizer("lass"))
    })

test_that("Distances in vectors",
    {
        a <- c(1, 1, 1, 1)
        b <- c(1, 0.8, 1, 1)
        expect_equal(mean_relative_tol(a, b), 0.05, tolerance = 0.001)

        expect_equal(mean_absolute_tol(a, b), 0.05, tolerance = 0.001)
    })

