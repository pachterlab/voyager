# Unit test validity function for SFEMethod class
library(spdep)

test_that("Wrong arguments in fun", {
    expect_error(SFEMethod(name = "sp.correlogram", variate = "uni",
                           scope = "global", package = "spdep",
                           title = "Correlogram", default_attr = NA,
                           fun = spdep::sp.correlogram,
                           reorganize_fun = Voyager:::.moran2df),
                 "The first two arguments of slot `fun` must be 'x' and 'listw'")
    expect_error(SFEMethod(name = "sp.correlogram", variate = "uni",
                           scope = "global", package = "spdep",
                           title = "Correlogram", default_attr = NA,
                           fun = spdep::sp.correlogram, use_graph = FALSE,
                           reorganize_fun = Voyager:::.moran2df),
                 "The first two arguments of slot `fun` must be 'x' and 'coords_df'")
})

test_that("Must have zero.policy", {
    expect_error(SFEMethod(name = "moran", variate = "uni",
                           scope = "global", package = "spdep",
                           title = "Moran's I", default_attr = NA,
                           fun = function(x, listw, ...) x,
                           reorganize_fun = Voyager:::.moran2df),
                 "zero.policy must be an argument of slot `fun`")
})

test_that("Check arguments of reorganize_fun", {
    expect_error(SFEMethod(name = "moran", variate = "uni",
                           scope = "global", package = "spdep",
                           title = "Moran's I", default_attr = NA,
                           fun = spdep::moran,
                           reorganize_fun = function(x) x),
                 "Slot `reorganize_fun` must have arguments 'out', 'name', and '...'")
    expect_error(SFEMethod(name = "localmoran", variate = "uni",
                           scope = "local", package = "spdep",
                           title = "Local Moran's I", default_attr = "Ii",
                           fun = spdep::localmoran,
                           reorganize_fun = function(x) x),
                 "Slot `reorganize_fun` must have arguments 'out', 'nb', and 'p.adjust.method'"
    )
})
