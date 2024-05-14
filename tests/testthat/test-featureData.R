library(SFEData)
library(SpatialFeatureExperiment)
library(scater)
library(spdep)

sfe <- McKellarMuscleData("small")
colGraph(sfe, "visium") <- findVisiumGraph(sfe)

# Already tested in test-univariate.R before I moved the functions to
# featureData.R, here just to test the edge case messages

test_that("Different packages", {
    sfe <- colDataMoransI(sfe, "nCounts")
    moran2 <- SFEMethod(
        name = "moran", title = "Moran's I", scope = "global",
        default_attr = NA,
        package = "scater", variate = "uni",
        fun = function(x, listw, zero.policy = NULL)
            spdep::moran(x, listw, n = length(listw$neighbours), S0 = Szero(listw),
                         zero.policy = zero.policy),
        reorganize_fun = Voyager:::.moran2df
    )
    expect_message(sfe <- colDataUnivariate(sfe, type = moran2,
                                            features = "nGenes", name = "moran"),
                   "please verify consistency between packages.")
})

test_that("Different versions of the same package", {
    sfe <- colDataMoransI(sfe, "nCounts")
    # Current version of spdep is 1.2.8
    metadata(colData(sfe))$params$moran$version <- package_version("1.2.6")
    expect_message(sfe <- colDataMoransI(sfe, "nGenes"),
                   "please verify consistency between versions.")
})

test_that("Same package and version but different parameters", {
    sfe <- colDataUnivariate(sfe, "sp.correlogram", features = "nCounts",
                             order = 3, style = "W")
    expect_error(sfe <- colDataUnivariate(sfe, "sp.correlogram",
                                          features = "nGenes", order = 3,
                                          style = "C"),
                 "please use a different name for different parameters.")
})

test_that("Different packages for graph", {
    sfe <- colDataMoransI(sfe, "nCounts", colGraphName = "visium",
                          name = "moran")
    g <- colGraph(sfe, "visium")
    attr(g, "method")$package[[1]] <- "foobar"
    colGraph(sfe, "visium2") <- g
    expect_message(sfe <- colDataMoransI(sfe, "nGenes", colGraphName = "visium2",
                                         name = "moran"),
                   "New results used package foobar")
})

test_that("Different versions of the same package for graph", {
    sfe <- colDataMoransI(sfe, "nCounts", colGraphName = "visium")
    g <- colGraph(sfe, "visium")
    attr(g, "method")$package[[2]] <- package_version("1.0.8")
    colGraph(sfe, "visium2") <- g
    expect_message(sfe <- colDataMoransI(sfe, "nGenes", colGraphName = "visium2"),
                   "New results used version 1.0.8")
})

test_that("Different functions used for graph", {
    sfe <- colDataMoransI(sfe, "nCounts", colGraphName = "visium")
    colGraph(sfe, "knn") <- findSpatialNeighbors(sfe, MARGIN = 2, method = "knearneigh",
                                                 k = 5, nn_method = "spdep")
    expect_message(sfe <- colDataMoransI(sfe, "nGenes", colGraphName = "knn"),
                   "New results used a spatial neighborhood graph computed with function")
})

test_that("Different parameters of the same function for graph", {
    colGraph(sfe, "knn") <- findSpatialNeighbors(sfe, MARGIN = 2, method = "knearneigh",
                                                 k = 5, nn_method = "spdep")
    sfe <- colDataMoransI(sfe, "nCounts", colGraphName = "knn")
    colGraph(sfe, "knn2") <- findSpatialNeighbors(sfe, MARGIN = 2, method = "knearneigh",
                                                  k = 6, nn_method = "spdep",
                                                  dist_type = "idw")
    expect_message(sfe <- colDataMoransI(sfe, "nGenes", colGraphName = "knn2"),
                   "was computed with different parameters")
})
