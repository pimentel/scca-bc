context("SCCAB interface")

test_that("sccab_subsample throws exception",
    {
        X <- matrix(0, nrow = 100, ncol = 10)
        expect_error(sccab_subsample(X, lam = 3, prop = 1))
        expect_error(sccab_subsample(X, lam = 3, prop = 0))
        expect_error(sccab_subsample(X, lam = 3, prop = 1.3))
        expect_error(sccab_subsample(X, lam = 3, prop = 0.001))
        expect_error(sccab_subsample(X, lam = 3, prop = -10))

        expect_error(sccab_subsample(data.frame(X), lam = 3, prop = -10))
    })

test_that("sccab throws exception",
    {
        X <- data.frame(matrix(0, nrow = 100, ncol = 10))
        expect_error(sccab(X, lam = 3))
    })
