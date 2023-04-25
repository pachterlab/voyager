library(SFEData)
library(scater)
library(scran)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(vdiffr)

sfe <- McKellarMuscleData("small")
sfe <- sfe[,sfe$in_tissue]
sfe <- logNormCounts(sfe)
gs <- modelGeneVar(sfe)
hvgs <- getTopHVGs(gs, fdr.threshold = 0.01)

g <- colGraph(sfe, "visium") <- findVisiumGraph(sfe)

mat_x <- logcounts(sfe)[hvgs[1:3],]
rownames(mat_x) <- rowData(sfe)[hvgs[1:3], "symbol"]
mat_y <- logcounts(sfe)[hvgs[4:7],]
rownames(mat_y) <- rowData(sfe)[hvgs[4:7], "symbol"]

test_that("Errors", {
    expect_error(calculateBivariate(sfe$nCounts, type = "lee", listw = g),
                 "y must be specified for vector x.")
    expect_error(calculateBivariate(sfe$nCounts, sfe$nGenes[-1], type = "lee",
                                    listw = g),
                 "same number")
    expect_error(calculateBivariate(mat_x, mat_y[,-1], type = "lee", listw = g),
                 "same number")
    mat_x2 <- mat_x
    rownames(mat_x2) <- NULL
    expect_error(calculateBivariate(mat_x2, mat_y, type = "lee.test", listw = g),
                 "Matrices x and y must have row names.")
})

test_that("Lee, taking matrix input", {
    out <- calculateBivariate(mat_x, type = "lee", listw = g)
    # Consistent with spdep
    ref <- spdep::lee(mat_x[1,], mat_x[2,], listw = g, n = ncol(mat_x))
    expect_equal(out[1,2], ref$L)
    # Another edge weight scheme
    g2 <- findVisiumGraph(sfe, style = "B")
    out2 <- calculateBivariate(mat_x, type = "lee", listw = g2)
    ref2 <- spdep::lee(mat_x[1,], mat_x[2,], listw = g2, n = ncol(mat_x))
    expect_equal(out2[1,2], ref2$L)

    out3 <- calculateBivariate(mat_x[1,], mat_y[1,], listw = g, type = "lee")
    ref3 <- spdep::lee(mat_x[1,], mat_y[1,], listw = g, n = ncol(mat_x))
    expect_equal(out3, ref3$L)

    out4 <- calculateBivariate(mat_x, mat_y, listw = g, type = "lee")
    expect_equal(rownames(out4), rownames(mat_x))
    expect_equal(colnames(out4), rownames(mat_y))
})

names_expect_mc <- c(
    "statistic", "parameter", "p.value", "alternative",
    "method", "res"
)
names_expect_mc <- paste("lee.mc", names_expect_mc, sep = "_")
test_that("lee.mc gives appropriate output", {
    # Vector x and vector y
    out <- calculateBivariate(mat_x[1,], mat_y[1,], listw = g, type = "lee.mc",
                              nsim = 49)
    expect_s4_class(out, "DataFrame")
    # Matrix x
    out2 <- calculateBivariate(mat_x, listw = g, type = "lee.mc", nsim = 49)
    expect_s4_class(out2, "DataFrame")
    expect_equal(names(out2), names_expect_mc)
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_x))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(rownames(out2), rns_expect)

    # Matrix x and matrix y
    out3 <- calculateBivariate(mat_x, mat_y, listw = g, type = "lee.mc", nsim = 49)
    expect_s4_class(out3, "DataFrame")
    expect_equal(names(out3), names_expect_mc)
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_y))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(rownames(out3), rns_expect)
})

names_expect_lt <- c(
    "statistic", "p.value", "alternative", "method",
    "Lee.s.L.statistic", "Expectation", "Variance"
)
names_expect_lt <- paste("lee.test", names_expect_lt, sep = "_")
test_that("lee.test gives appropriate results", {
    out <- calculateBivariate(mat_x[1,], mat_y[1,], listw = g, type = "lee.test")
    expect_s4_class(out, "DataFrame")

    out2 <- calculateBivariate(mat_x, mat_y, listw = g, type = "lee.test")
    expect_s4_class(out2, "DataFrame")
    expect_equal(names(out2), names_expect_lt)
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_y))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(rownames(out2), rns_expect)
})

test_that("Local Lee results", {
    out <- calculateBivariate(mat_x, mat_y, listw = g, type = "locallee")
    expect_s4_class(out, "DataFrame")
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_y))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(names(out), rns_expect)
})

names_expect_lb <- c("Ibvi", "E.Ibvi", "Var.Ibvi", "Z.Ibvi", "Pr(z != E(Ibvi))",
                     "Pr(z != E(Ibvi)) Sim", "Pr(folded) Sim", "-log10p Sim",
                     "-log10p_adj Sim")
test_that("localmoral_bv results", {
    out <- calculateBivariate(mat_x, mat_y, listw = g, type = "localmoran_bv",
                              nsim = 49)
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_y))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(names(out), rns_expect)
    expect_true(all(vapply(out, is.matrix, FUN.VALUE = logical(1))))
    expect_true(all(vapply(out, function(o) identical(colnames(o), names_expect_lb),
                           FUN.VALUE = logical(1))))
})

df <- df2sf(spatialCoords(sfe), spatialCoordsNames(sfe))
test_that("Cross variogram output", {
    expect_message(out <- calculateBivariate(mat_x, mat_y, type = "cross_variogram",
                                             coords_df = df),
                   "Cross correlograms within columns of x and within columns of y are also computed.")
    out <- calculateBivariate(mat_x, type = "cross_variogram", coords_df = df)
    expect_s3_class(out, "data.frame")
    expect_named(out, c("np", "dist", "gamma", "dir.hor", "dir.ver", "id"))
    expect_doppelganger("Plot cross variograms", {
        plotCrossVariogram(out)
    })
    expect_doppelganger("Plot cross variograms no np", {
        plotCrossVariogram(out, show_np = FALSE)
    })
    # Multiple angles
    out2 <- calculateBivariate(mat_x, type = "cross_variogram", coords_df = df,
                               alpha = c(30, 90, 150))
    expect_doppelganger("Plot cross variograms with anisotropy", {
        plotCrossVariogram(out2)
    })
})

test_that("Cross variogram map", {
    out <- calculateBivariate(mat_x, type = "cross_variogram_map", width = 500,
                              cutoff = 2000, coords_df = df)
    expect_s3_class(out, "data.frame")
    expect_named(out, c("Myh2.Myh1", "np.Myh2.Myh1", "Csrp3.Myh1", "np.Csrp3.Myh1",
                        "Myh1", "np.Myh1", "Myh2.Csrp3", "np.Myh2.Csrp3", "Csrp3", "np.Csrp3",
                        "Myh2", "np.Myh2", "dx", "dy"))
    expect_doppelganger("Plot cross variogram map", {
        plotCrossVariogramMap(out)
    })
})

test_that("Cross variogram map for one pair", {
    out <- calculateBivariate(x = mat_x[1,], y = mat_y[1,],
                              type = "cross_variogram_map", width = 500,
                              cutoff = 2000, coords_df = df)
    expect_named(out, c("x.y", "np.x.y", "y", "np.y", "x", "np.x", "dx", "dy"))
})

test_that("calculateBivariate SFE method", {
    expect_error(calculateBivariate(sfe, type = "lee", feature1 = hvgs[1]),
                 "feature2 must be specified when feature1 has length 1.")
    out <- calculateBivariate(sfe, type = "lee.test", feature1 = hvgs[1:2],
                              feature2 = hvgs[3:4], colGraphName = "visium")
    expect_s4_class(out, "DataFrame")
    expect_equal(names(out), names_expect_lt)
    rns_expect <- expand.grid(hvgs[1:2], hvgs[3:4])
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(rownames(out), rns_expect)

    out2 <- calculateBivariate(sfe, type = "cross_variogram", feature1 = hvgs[1:2],
                               colGeometryName = "spotPoly")
    expect_s3_class(out2, "data.frame")
    expect_named(out2, c("np", "dist", "gamma", "dir.hor", "dir.ver", "id"))
})

sfe2 <- McKellarMuscleData("small2")
sfe2 <- sfe2[,sfe2$in_tissue]
sfe2 <- logNormCounts(sfe2)
sfe3 <- SpatialFeatureExperiment::cbind(sfe, sfe2)
colGraphs(sfe3, sample_id = "all", name = "visium") <- findVisiumGraph(sfe3, sample_id = "all")

test_that("calculateBivariate SFE method, two samples", {
    out <- calculateBivariate(sfe3, type = "lee.test", feature1 = hvgs[1:2],
                              feature2 = hvgs[3:4], colGraphName = "visium",
                              sample_id = "all")
    expect_type(out, "list")
    expect_named(out, sampleIDs(sfe3))
    expect_s4_class(out[[1]], "DataFrame")
    expect_equal(names(out[[1]]), names_expect_lt)
    rns_expect <- expand.grid(hvgs[1:2], hvgs[3:4])
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(rownames(out[[1]]), rns_expect)

    out2 <- calculateBivariate(sfe3, type = "cross_variogram", feature1 = hvgs[1:2],
                               colGeometryName = "spotPoly", sample_id = "all")
    expect_type(out2, "list")
    expect_named(out2, sampleIDs(sfe3))
    expect_s3_class(out2[[1]], "data.frame")
    expect_named(out2[[1]], c("np", "dist", "gamma", "dir.hor", "dir.ver", "id"))
})

test_that("runBivariate", {
    expect_error(runBivariate(sfe, "lee", feature1 = rownames(mat_x)),
                 "Global bivariate results can't be stored in the SFE object.")
    sfe <- runBivariate(sfe, "locallee", feature1 = rownames(mat_x),
                        feature2 = rownames(mat_y), colGraphName = "visium",
                        swap_rownames = "symbol")
    rns_expect <- expand.grid(rownames(mat_x), rownames(mat_y))
    rns_expect <- paste(rns_expect[[1]], rns_expect[[2]], sep = "__")
    expect_equal(localResultFeatures(sfe, "locallee"), rns_expect)
})
