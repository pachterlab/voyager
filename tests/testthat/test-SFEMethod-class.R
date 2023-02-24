# Unit test validity function for SFEMethod class
library(spdep)

test_that("Wrong info names", {
    expect_error(SFEMethod(info = c(foo = "bar"), fun = spdep::moran,
                           reorganize_fun = Voyager:::.moran2df),
                 "Slot `info` must have names")
})

test_that("Wrong variates", {
    expect_error(SFEMethod(info = c(name = "moran", variate = "foo",
                                    scope = "global", package = "spdep",
                                    title = "Moran's I", default_attr = NA),
                           fun = spdep::moran,
                           reorganize_fun = Voyager:::.moran2df),
                 "Field 'variate' in slot `info` must be one of")
})

test_that("Wrong scope", {
    expect_error(SFEMethod(info = c(name = "moran", variate = "uni",
                                    scope = "bar", package = "spdep",
                                    title = "Moran's I", default_attr = NA),
                           fun = spdep::moran,
                           reorganize_fun = Voyager:::.moran2df),
                 "Field 'scope' in slot `info` must be one of")
})

test_that("Wrong arguments in fun", {
    expect_error(SFEMethod(info = c(name = "sp.correlogram", variate = "uni",
                                    scope = "global", package = "spdep",
                                    title = "Correlogram", default_attr = NA),
                           fun = spdep::sp.correlogram,
                           reorganize_fun = Voyager:::.moran2df),
                 "The first two arguments of slot `fun` must be 'x' and 'listw'")
})

test_that("Must have zero.policy", {
    expect_error(SFEMethod(info = c(name = "moran", variate = "uni",
                                    scope = "global", package = "spdep",
                                    title = "Moran's I", default_attr = NA),
                           fun = function(x, listw, ...) x,
                           reorganize_fun = Voyager:::.moran2df),
                 "zero.policy must be an argument of slot `fun`")
})

test_that("Check arguments of reorganize_fun", {
    expect_error(SFEMethod(info = c(name = "moran", variate = "uni",
                                    scope = "global", package = "spdep",
                                    title = "Moran's I", default_attr = NA),
                           fun = spdep::moran,
                           reorganize_fun = function(x) x),
                 "Slot `reorganize_fun` must have arguments 'out', 'name', and '...'")
    expect_error(SFEMethod(info = c(name = "localmoran", variate = "uni",
                                    scope = "local", package = "spdep",
                                    title = "Local Moran's I", default_attr = "Ii"),
                           fun = spdep::localmoran,
                           reorganize_fun = function(x) x),
                 "Slot `reorganize_fun` must have arguments 'out', 'nb', and 'p.adjust.method'"
    )
})
