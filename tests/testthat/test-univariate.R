library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(Matrix)
library(bluster)
library(scater)
sfe <- readRDS(system.file("extdata/sfe.rds", package = "Voyager"))
set.seed(29)
mat <- assay(sfe, "counts")
mat1 <- mat[, colData(sfe)$sample_id == "sample01"]

out_m <- calculateMoransI(mat1, listw = colGraph(sfe, "visium", sample_id = "sample01"))
test_that("Error when a multivariate method is used", {
    expect_error({
        calculateUnivariate(mat1, type = "multispati",
                            listw = colGraph(sfe, "visium", sample_id = "sample01"))
    }, "`type` must be a univariate method.")
})

test_that("Correct structure of calculateMoransI output (matrix)", {
    expect_s4_class(out_m, "DataFrame")
    expect_equal(names(out_m), c("moran", "K"))
    expect_true(is.numeric(out_m$moran))
    expect_true(is.numeric(out_m$K))
    expect_equal(rownames(out_m), rownames(mat1))
})

test_that("Correctly add results to rowData when features = NULL", {
    sfe <- runUnivariate(sfe, type = "geary", features = NULL,
                         colGraphName = "visium", sample_id = "sample01",
                         exprs_values = "counts")
    expect_true(is.numeric(rowData(sfe)$geary_sample01))
    expect_true(all(!is.na(rowData(sfe)$geary_sample01)))
})

test_that("Correctly add results to rowData when features = NULL with Moran's I wrapper", {
    sfe <- runMoransI(sfe, features = NULL, colGraphName = "visium",
                      sample_id = "sample01", exprs_values = "counts")
    expect_true(is.numeric(rowData(sfe)$moran_sample01))
    expect_true(all(!is.na(rowData(sfe)$moran_sample01)))
})

test_that("Correct structure of colDataMoransI output", {
    out <- colDataMoransI(sfe, "nCounts", "visium", sample_id = "sample01")
    expect_s4_class(out, "SpatialFeatureExperiment")
    fd <- colFeatureData(out)
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), c("moran_sample01", "K_sample01"))
    expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
    expect_true(is.na(fd["barcode", "moran_sample01"]))
    expect_false(is.na(fd["nCounts", "moran_sample01"]))
    # Check the params field
    params <- getParams(out, "moran", colData = TRUE)
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$graph_params,
                 attr(colGraph(out, "visium", "sample01"), "method"))
})

test_that("Correct structure of colGeometryMoransI output", {
    out <- colGeometryMoransI(sfe,
        colGeometryName = "spotPoly",
        colGraphName = "visium", features = "foo",
        sample_id = "sample01"
    )
    fd <- geometryFeatureData(out, "spotPoly")
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), c("moran_sample01", "K_sample01"))
    expect_equal(rownames(fd), c("geometry", "foo"))
    expect_true(is.na(fd["geometry", "moran_sample01"]))
    expect_false(is.na(fd["foo", "moran_sample01"]))
    # Check the params field
    params <- getParams(out, "moran", colGeometryName = "spotPoly")
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$graph_params,
                 attr(colGraph(out, "visium", "sample01"), "method"))
})

test_that("Correct structure of calculateMoransI output (SFE)", {
    out <- calculateMoransI(sfe,
        colGraphName = "visium", features = rownames(mat1),
        sample_id = "sample01", exprs_values = "counts"
    )
    expect_s4_class(out, "DataFrame")
    expect_equal(names(out), c("moran", "K"))
    expect_true(is.numeric(out$moran))
    expect_true(is.numeric(out$K))
    expect_equal(rownames(out), rownames(mat1))
})

test_that("Properly add Moran's I results (no permutation) to SFE rowData", {
    sfe2 <- runMoransI(sfe,
        colGraphName = "visium", features = rownames(mat1),
        sample_id = "sample01",
        exprs_values = "counts"
    )
    rd <- rowData(sfe2)
    expect_equal(names(rd), c("moran_sample01", "K_sample01"))
    names(rd) <- c("moran", "K")
    # just check the values, checking attributes a little later
    metadata(rd) <- list()
    expect_equal(rd, out_m)
    # Check the params field
    params <- getParams(sfe2, "moran")
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$graph_params,
                 attr(colGraph(sfe2, "visium", "sample01"), "method"))
})

names_expect_mc <- c(
    "statistic", "parameter", "p.value", "alternative",
    "method", "data.name", "res"
)
names_expect_mc <- paste("moran.mc", names_expect_mc, sep = "_")
names_expect_mc_sample <- paste(names_expect_mc, "sample01", sep = "_")

test_that("Correct structure of mc.sim output", {
    out <- calculateUnivariate(mat1,
        listw = colGraph(sfe, "visium", sample_id = "sample01"),
        type = "moran.mc", nsim = 10
    )
    expect_s4_class(out, "DataFrame")
    expect_equal(names(out), names_expect_mc)
    expect_equal(rownames(out), rownames(mat1))
    expect_equal(unname(out$moran.mc_statistic), out_m$moran)
})

test_that("Correct structure of colDataUnivariate output, with list column", {
    out <- colDataUnivariate(sfe,
        type = "moran.mc", colGraphName = "visium",
        features = "nCounts", sample_id = "sample01",
        nsim = 10
    )
    fd <- colFeatureData(out)
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), names_expect_mc_sample)
    expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
    expect_true(is.na(fd["barcode", "moran.mc_statistic_sample01"]))
    expect_false(is.na(fd["nCounts", "moran.mc_statistic_sample01"]))
})

test_that("Correct structure of colGeometryUnivariate output, with list column", {
    out <- colGeometryUnivariate(sfe,
        colGeometryName = "spotPoly",
        type = "moran.mc",
        colGraphName = "visium", features = "foo",
        sample_id = "sample01", nsim = 10
    )
    fd <- geometryFeatureData(out, "spotPoly")
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), names_expect_mc_sample)
    expect_equal(rownames(fd), c("geometry", "foo"))
    expect_true(is.na(fd["geometry", "moran.mc_statistic_sample01"]))
    expect_false(is.na(fd["foo", "moran.mc_statistic_sample01"]))
})

test_that("MoranMC results properly added to rowData", {
    sfe2 <- runUnivariate(sfe,
        type = "moran.mc", colGraphName = "visium",
        features = rownames(mat1),
        sample_id = "sample01", exprs_values = "counts",
        nsim = 10
    )
    rd <- rowData(sfe2)
    expect_equal(names(rd), names_expect_mc_sample)
    expect_equal(rd$moran.mc_statistic_sample01, out_m$moran)
})

test_that("DataFrame results for sp.correlogram, method = I", {
    out <- calculateUnivariate(mat1,
        type = "sp.correlogram",
        listw = colGraph(sfe, "visium",
            sample_id = "sample01"
        ),
        order = 2, method = "I"
    )
    expect_s4_class(out, "DFrame")
    expect_equal(rownames(out), rownames(mat1))
    i1 <- vapply(out[, 1], function(o) o[1, 1], FUN.VALUE = numeric(1))
    expect_equal(unname(i1), out_m$moran)
})

test_that("DataFrame results for sp.correlogram, method = corr", {
    out <- calculateUnivariate(mat1,
                               type = "sp.correlogram",
                               listw = colGraph(sfe, "visium",
                                                sample_id = "sample01"
                               ),
                               order = 2, method = "corr"
    )
    expect_s4_class(out, "DFrame")
    expect_equal(rownames(out), rownames(mat1))
    expect_true(all(vapply(out[,1], function(x) is.numeric(x) & length(x) == 2L,
                           FUN.VALUE = logical(1))))
})

names_expect_mp <- c(
    "x", "wx", "is_inf", "labels", "dfb.1_", "dfb.x",
    "dffit", "cov.r", "cook.d", "hat"
)

test_that("DataFrame output for moran.plot, and some local results in general", {
    out <- calculateUnivariate(mat1,
        listw = colGraph(sfe, "visium", sample_id = "sample01"),
        type = "moran.plot"
    )
    expect_s4_class(out, "DFrame")
    expect_true(all(vapply(out, is.data.frame, FUN.VALUE = logical(1))))
    expect_equal(names(out[[1]]), names_expect_mp)
    expect_equal(names(out), rownames(mat1))
    expect_equal(nrow(out), ncol(mat1))
})

names_expect_gg <- c(
    "statistic", "p.value", "alternative", "data.name", "method",
    "Global.G.statistic", "Expectation", "Variance"
)
names_expect_gg <- paste("globalG.test", names_expect_gg, sep = "_")
test_that("DataFrame output for globalG.test, or htest in general", {
    listw <- colGraph(sfe, "visium", sample_id = "sample01")
    listw2 <- nb2listw(listw$neighbours, style = "B")
    out <- calculateUnivariate(mat1, listw = listw2, type = "globalG.test")
    expect_s4_class(out, "DFrame")
    expect_equal(names(out), names_expect_gg)
    expect_equal(rownames(out), rownames(mat1))
})

names_expect_lm <- c(
    "Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr(z != E(Ii))", "mean",
    "median", "pysal", "-log10p", "-log10p_adj"
)
test_that("DataFrame output for localmoran", {
    out <- calculateUnivariate(mat1,
        listw = colGraph(sfe, "visium", sample_id = "sample01"),
        type = "localmoran"
    )
    expect_s4_class(out, "DFrame")
    expect_true(all(vapply(out, is.data.frame, FUN.VALUE = logical(1))))
    expect_equal(names(out[[1]]), names_expect_lm)
    expect_equal(names(out), rownames(mat1))
    expect_equal(nrow(out), ncol(mat1))
})

names_expect_lg <- c(
    "localG", "Gi", "E.Gi", "Var.Gi", "StdDev.Gi", "Pr(z != E(Gi))",
    "Pr(z != E(Gi)) Sim", "Pr(folded) Sim", "Skewness",
    "Kurtosis", "-log10p Sim", "-log10p_adj Sim", "cluster"
)
test_that("DataFrame output for localG_perm", {
    out <- calculateUnivariate(mat1,
        listw = colGraph(sfe, "visium", sample_id = "sample01"),
        type = "localG_perm"
    )
    expect_s4_class(out, "DFrame")
    expect_true(all(vapply(out, is.data.frame, FUN.VALUE = logical(1))))
    expect_equal(colnames(out[[1]]), names_expect_lg)
    expect_equal(names(out), rownames(mat1))
    expect_equal(nrow(out), ncol(mat1))
})

#test_that("DataFrame output for localG, not perm, when output is a vector", {
#    out <- calculateUnivariate(mat1,
#                               listw = colGraph(sfe, "visium", "sample01"),
#                               type = "localG", return_internals = FALSE)
#    expect_s4_class(out, "DFrame")
#    expect_true(all(vapply(out, function(o) is.atomic(o) & is.vector(o) &
#                               is.numeric(o), FUN.VALUE = logical(1))))
#    expect_equal(names(out), rownames(mat1))
#    expect_equal(nrow(out), ncol(mat1))
#})

test_that("Properly add localmoran results to localResults when there're multiple samples", {
    feats_use <- rownames(mat1)[1:2]
    sfe2 <- runUnivariate(sfe,
                          type = "localmoran", colGraphName = "visium",
                          features = feats_use,
                          sample_id = "all", exprs_values = "counts"
    )
    expect_equal(localResultNames(sfe2), "localmoran")
    lrs <- localResults(sfe2, "localmoran", feats_use, sample_id = "all")
    expect_equal(names(lrs), feats_use)
    expect_true(all(!is.na(lrs[[1]]$Ii)))
    expect_true(is.numeric(lrs[[1]]$Ii))
    expect_true(all(!is.na(lrs[[2]]$Ii)))
    expect_true(is.numeric(lrs[[2]]$Ii))
})

library(SFEData)
sfe1 <- McKellarMuscleData("small")
sfe2 <- McKellarMuscleData("small2")
sfe <- SpatialFeatureExperiment::cbind(sfe1, sfe2)
colGraphs(sfe, name = "visium", sample_id = "all") <- findVisiumGraph(sfe, sample_id = "all")
sfe <- colDataUnivariate(sfe, "localmoran", features = "nCounts", sample_id = "all")
res <- localResult(sfe, "localmoran", "nCounts", sample_id = "all")

test_that("Properly add moran.plot results to localResults when only one gene is used", {
    colGraph(sfe1, "visium") <- findVisiumGraph(sfe1)
    sfe1 <- logNormCounts(sfe1)
    sfe1 <- runUnivariate(sfe1, "moran.plot", features = "Myh1", colGraphName = "visium",
                          swap_rownames = "symbol")
    expect_equal(localResultFeatures(sfe1, "moran.plot"),
                 .symbol2id(sfe1, "Myh1", "symbol"))
    lr <- localResult(sfe1, "moran.plot", "Myh1", swap_rownames = "symbol")
    expect_s3_class(lr, "data.frame")
    expect_equal(names(lr), names_expect_mp)
})

test_that("colDataUnivariate run on multiple samples", {
    expect_s3_class(res, "data.frame")
    expect_equal(names(res), names_expect_lm)
    # parameters
    params <- getParams(sfe, "localmoran", local = TRUE, colData = TRUE)
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$p.adjust.method, "BH")
    expect_equal(params$graph_params,
                 attr(colGraph(sfe, "visium", "Vis5A"), "method"))
})

test_that("colGeometryUnivariate run on multiple samples", {
    colGeometry(sfe, "spotPoly", sample_id = "all")$foo <- sfe$nCounts
    sfe <- colGeometryUnivariate(sfe, "localmoran", features = "foo",
                                 sample_id = "all", colGeometryName = "spotPoly")
    res2 <- localResult(sfe, "localmoran", "foo", colGeometryName = "spotPoly",
                        sample_id = "all")
    class(res2) <- "data.frame"
    expect_equal(res, res2, ignore_attr = "row.names")
})

annotGeometry(sfe, "myofiber_simplified", "all") <-
    sf::st_buffer(annotGeometry(sfe, "myofiber_simplified", "all"), dist = 0)


test_that("annotGeometryUnivariate run on multiple samples", {
    # No message about graph parameters.
    # Message comes from not getting params from annotGeometries.
    expect_message(annotGraphs(sfe, name = "knn", sample_id = "all") <-
                       findSpatialNeighbors(sfe, MARGIN = 3, type = "myofiber_simplified",
                                            method = "knearneigh", k = 5, sample_id = "all"),
                   regexp = NA)
    sfe <- annotGeometryUnivariate(sfe, "localmoran", features = "area",
                                   sample_id = "all",
                                   annotGeometryName = "myofiber_simplified")
    res3 <- localResult(sfe, "localmoran", "area",
                        annotGeometryName = "myofiber_simplified",
                        sample_id = "all")
    expect_s3_class(res, "data.frame")
    expect_equal(names(res), names_expect_lm)
    # parameters
    params <- getParams(sfe, "localmoran", local = TRUE,
                        annotGeometryName = "myofiber_simplified")
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$p.adjust.method, "BH")
    expect_equal(params$graph_params,
                 attr(annotGraph(sfe, "knn", "Vis5A"), "method"))
})

sfe <- logNormCounts(sfe)
sfe <- runPCA(sfe, ncomponent = 2)

test_that("Univariate global results corrected added to metadata of reducedDim", {
    sfe <- reducedDimMoransI(sfe, "PCA", components = 1:2, sample_id = "Vis5A")
    fd <- reducedDimFeatureData(sfe, "PCA")
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), c("moran_Vis5A", "K_Vis5A"))
    expect_equal(rownames(fd), c("PC1", "PC2"))
    # parameters
    params <- getParams(sfe, "moran", reducedDimName = "PCA")
    expect_equal(params$package, "spdep")
    expect_equal(params$version, packageVersion("spdep"))
    expect_null(params$zero.policy)
    expect_false(params$include_self)
    expect_equal(params$graph_params,
                 attr(colGraph(sfe, "visium", "Vis5A"), "method"))
})

test_that("Univariate local results for reducedDim", {
    sfe <- reducedDimUnivariate(sfe, "localmoran", dimred = "PCA", components = 1)
    expect_true("PC1" %in% localResultFeatures(sfe, "localmoran"))
    lr <- localResult(sfe, "localmoran", feature = "PC1", sample_id = "Vis5A")
    expect_equal(colnames(lr), names_expect_lm)
    expect_true(all(vapply(lr, is.numeric, FUN.VALUE = logical(1))))
})

test_that("Univariate global results corrected added after already running for something else", {
    sfe <- colDataUnivariate(sfe, "sp.correlogram", "nCounts", order = 5, zero.policy = TRUE)
    sfe <- reducedDimUnivariate(sfe, "sp.correlogram", order = 3)
    fd <- reducedDimFeatureData(sfe, "PCA")
    expect_true("sp.correlogram_I_Vis5A" %in% names(fd))
})

# to do: gstat, univariate method not using listw
