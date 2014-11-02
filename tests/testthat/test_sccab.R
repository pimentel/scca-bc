context("sccab_param")
test_that("valid creation of sccab params",
    {
        p <- sccab_params(10)

        expect_is(p, "sccab_params")
        expect_is(p$apply_fun, "function")
        expect_is(p$optim_fun, "function")

        expect_equal(p$d_upr, 10)
        expect_equal(p$optim_fun, lasso_max_d)

        expect_equal(sort(names(p)), sort(c("d_lwr", "d_upr", "ab_lam",
                    "n_samp", "prop", "apply_fun", "verbose", "optim_fun", "optim_str")))

        expect_equal(sccab_params(4)$d_upr, 4)

        expect_error(sccab_params(4, ab_lam = c(-1, 1, 3)))
        expect_error(sccab_params(4, ab_lam = c(0, 1, 3)))

        expect_error(sccab_params(4, optim_fun = 0.9))
        expect_error(sccab_params(4, optim_fun = "meow"))
        expect_equal(sccab_params(4, optim_fun = "timeseries")$optim_fun,
            timeseries_max_d)

        expect_error(sccab_params(3, prop = -1.0))

        expect_warning(sccab_params(3, prop = 1.0))

        expect_error(sccab_params(d_lwr = 4, d_upr = 3))
        expect_error(sccab_params(4, n_samp = 0.9))

        expect_error(sccab_params(4, n_samp = 0.9))
    })

context("sccab_result")
test_that("valid sccab result",
    {
        dummy_sccab <- lapply(1:50, function(x)
            list(ab = rnorm(100), d = rnorm(20),
                sccaLam = matrix(c(1,1), ncol = 2)))

        p <- sccab_params(20, n_samp = 50)
        res <- sccab_result(dummy_sccab, p)
        expect_is(res, "sccab_result")
    })
