library(SFEData)
library(scater)
library(SpatialExperiment)
library(vdiffr)
library(bluster)
library(Matrix)

expect_ggplot <- function(g) {
    expect_s3_class(g, "ggplot")
    expect_error(ggplot_build(g), NA)
}

sfe <- McKellarMuscleData("small")
sfe <- sfe[,sfe$in_tissue]
sfe <- logNormCounts(sfe)
gs <- order(Matrix::rowSums(counts(sfe)), decreasing = TRUE)[1:10]
genes <- rownames(sfe)[gs]
mat <- logcounts(sfe)[gs,]
df <- df2sf(spatialCoords(sfe), spatialCoordsNames(sfe))

test_that("calculateUnivariate (matrix) returns correct output for variogram", {
    out <- calculateUnivariate(mat, variogram, coords_df = df, model = "Sph")
    expect_s4_class(out, "DataFrame")
    expect_equal(nrow(out), 10)
    expect_named(out, "variogram")
    expect_type(out$variogram, "list")
    expect_true(all(vapply(out$variogram, function(x) is(x, "autofitVariogram"),
                           FUN.VALUE = logical(1))))
})

test_that("calculateUnivariate (matrix) returns correct output for variogram map", {
    out <- calculateUnivariate(mat, variogram_map, coords_df = df,
                               cutoff = 2000, width = 500)
    expect_s4_class(out, "DataFrame")
    expect_equal(nrow(out), 10)
    expect_named(out, "variogram_map")
    expect_type(out$variogram_map, "list")
    expect_true(all(vapply(out$variogram, function(x) is.data.frame(x),
                           FUN.VALUE = logical(1))))
    dfs <- do.call(rbind, out$variogram_map)
    expect_named(dfs, c("var1", "np.var1", "dx", "dy"))
})

test_that("calculateUnivariate (SFE) returns correct output for variogram", {
    out <- calculateUnivariate(sfe, variogram, features = rownames(sfe)[gs],
                               model = "Sph")
    expect_s4_class(out, "DataFrame")
    expect_equal(nrow(out), 10)
    expect_named(out, "variogram")
    expect_type(out$variogram, "list")
    expect_true(all(vapply(out$variogram, function(x) is(x, "autofitVariogram"),
                           FUN.VALUE = logical(1))))
})

test_that("colDataUnivariate for variogram (or when not using graph)", {
    sfe <- colDataUnivariate(sfe, variogram, features = c("nCounts", "nGenes"),
                             model = "Sph")
    out <- colFeatureData(sfe)
    expect_named(out, "variogram_Vis5A")
    expect_type(out$variogram_Vis5A, "list")
    expect_true(all(vapply(out[c("nCounts", "nGenes"), "variogram_Vis5A"],
                           function(x) is(x, "autofitVariogram"),
                           FUN.VALUE = logical(1))))
    expect_true(all(vapply(out[setdiff(names(colData(sfe)), c("nCounts", "nGenes")),
                               "variogram_Vis5A"],
                           is.na, FUN.VALUE = logical(1))))
    # Check old params when not using listw
    sfe <- colDataUnivariate(sfe, variogram, features = "prop_mito",
                             model = "Sph")
    out <- colFeatureData(sfe)
    expect_true(all(vapply(out[c("nCounts", "nGenes", "prop_mito"), "variogram_Vis5A"],
                           function(x) is(x, "autofitVariogram"),
                           FUN.VALUE = logical(1))))
    expect_true(all(vapply(out[setdiff(names(colData(sfe)),
                                       c("nCounts", "nGenes", "prop_mito")),
                               "variogram_Vis5A"],
                           is.na, FUN.VALUE = logical(1))))
    expect_message(sfe <- colDataUnivariate(sfe, variogram, features = c("nCounts", "nGenes"),
                                            model = "Sph", alpha = c(30, 90, 150),
                                            name = "variogram_anis"),
                   "gstat does not fit anisotropic variograms")
})

sfe2 <- McKellarMuscleData("small2")
sfe2 <- sfe2[,sfe2$in_tissue]
sfe2 <- logNormCounts(sfe2)
sfe3 <- SpatialFeatureExperiment::cbind(sfe, sfe2)
sfe <- runPCA(sfe, ncomponents = 5)

test_that("plotVariogram", {
    sfe <- colDataUnivariate(sfe, variogram, features = c("nCounts", "nGenes"),
                             model = "Sph")
    expect_ggplot({
        plotVariogram(sfe, "nCounts")
    })
    expect_ggplot({
        plotVariogram(sfe, "nCounts", group = "features")
    })
    sfe <- reducedDimUnivariate(sfe, variogram, components = 1:3, model = "Ste")
    expect_ggplot({
        plotVariogram(sfe, 1:3, group = "features", reducedDimName = "PCA")
    })
    expect_ggplot({
        plotVariogram(sfe, 3:5, group = "features", reducedDimName = "PCA")
    })
    sfe3 <- colDataUnivariate(sfe3, variogram, features = c("nCounts", "nGenes"),
                              model = "Sph")
    expect_ggplot({
        plotVariogram(sfe3, "nCounts")
    })
    expect_ggplot({
        plotVariogram(sfe3, c("nCounts", "nGenes"), group = "feature")
    })
    expect_ggplot({
        plotVariogram(sfe3, c("nCounts", "nGenes"), group = "sample_id")
    })
    expect_ggplot({
        plotVariogram(sfe3, c("nCounts", "nGenes"))
    })
    color_values <- c("A", "B")
    expect_ggplot( {
        plotVariogram(sfe3, c("nCounts", "nGenes"), color_by = color_values)
    })
    expect_ggplot({
        plotVariogram(sfe3, c("nCounts", "nGenes"), color_by = color_values,
                      group = "features")
    })
    sfe <- colDataUnivariate(sfe, variogram, features = c("nCounts", "nGenes"),
                             model = "Sph", alpha = c(30, 90, 150),
                             name = "variogram_anis")
    expect_ggplot({
        plotVariogram(sfe, "nGenes", name = "variogram_anis")
    })
    expect_ggplot({
        plotVariogram(sfe, "nGenes", group = "angles", name = "variogram_anis")
    })
})

sfe <- runUnivariate(sfe, variogram, features = genes)
sfe3 <- runUnivariate(sfe3, variogram, features = genes)
test_that("clusterVariogram", {
    var_clusts <- clusterVariograms(sfe, features = genes,
                                    BLUSPARAM = HclustParam(),
                                    swap_rownames = "symbol")
    expect_s3_class(var_clusts, "data.frame")
    expect_named(var_clusts, c("feature", "cluster", "sample_id"))
    expect_s3_class(var_clusts$cluster, "factor")
    # Plotting clusters
    expect_ggplot({
        plotVariogram(sfe, genes[1:5], color_by = var_clusts, group = "features",
                      swap_rownames = "symbol")
    })
    expect_ggplot({
        plotVariogram(sfe, genes, color_by = var_clusts, group = "features",
                      use_lty = FALSE, swap_rownames = "symbol", show_np = FALSE)
    })
    # Across 2 samples
    var_clusts2 <- clusterVariograms(sfe3, features = genes,
                                     BLUSPARAM = HclustParam(),
                                     swap_rownames = "symbol")
    expect_ggplot({
        plotVariogram(sfe3, genes, color_by = var_clusts2, group = "features",
                      use_lty = FALSE, swap_rownames = "symbol", show_np = FALSE)
    })
    # When there's only one cluster, don't use color_by
    var_clusts3 <- var_clusts
    var_clusts3$cluster <- 1
    expect_ggplot({
        plotVariogram(sfe, genes[1:5], color_by = var_clusts3, group = "features",
                      swap_rownames = "symbol")
    })
    expect_message(plotVariogram(sfe, genes[1:5], color_by = var_clusts3, group = "features",
                                 swap_rownames = "symbol"),
                   "Only one cluster is present. Not using color_by.")
})

test_that("plotVariogramMap", {
    sfe3 <- colDataUnivariate(sfe3, variogram_map, c("nCounts", "nGenes"),
                              width = 500, cutoff = 2000)
    expect_doppelganger("One sample, one feature, variogram map", {
        plotVariogramMap(sfe3, "nCounts", sample_id = "Vis5A")
    })
    expect_doppelganger("Plot np", {
        plotVariogramMap(sfe3, "nCounts", sample_id = "Vis5A", plot_np = TRUE)
    })
    expect_doppelganger("Multiple samples, one feature", {
        plotVariogramMap(sfe3, "nCounts")
    })
    expect_doppelganger("Multiple samples, multiple features", {
        plotVariogramMap(sfe3, c("nCounts", "nGenes"))
    })
})
