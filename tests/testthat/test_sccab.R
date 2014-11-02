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

context("SCCAB parameters")
test_that("valid creation of sccab params",
    {
        p <- sccab_params(10)

        expect_is(p, "sccab_params")
        expect_is(p$apply_fun, "function")
        expect_is(p$optim_fun, "function")

        expect_equal(p$d_upr, 10)
        expect_equal(p$optim_fun, lasso_max_d)

        expect_equal(sort(names(p)), sort(c("d_lwr", "d_upr", "ab_lam",
                    "n_samp", "prop", "apply_fun", "verbose", "optim_fun")))

        expect_equal(sccab_params(4)$d_upr, 4)

        expect_error(sccab_params(4, ab_lam = c(-1, 1, 3)))
        expect_error(sccab_params(4, ab_lam = c(0, 1, 3)))

        expect_error(sccab_params(4, optim_fun = 0.9))
        expect_error(sccab_params(4, optim_fun = "meow"))
        expect_equal(sccab_params(4, optim_fun = "timeseries")$optim_fun,
            timeseries_max_d)

        expect_warning(sccab_params(3, prop = 1.0))

        expect_error(sccab_params(d_lwr = 4, d_upr = 3))
        expect_error(sccab_params(4, n_samp = 0.9))

        expect_error(sccab_params(4, n_samp = 0.9))
    })
