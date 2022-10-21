library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(Matrix)
library(bluster)
sfe <- readRDS(system.file("extdata/sfe.rds", package = "Voyager"))
set.seed(29)
mat <- assay(sfe, "counts")
mat1 <- mat[, colData(sfe)$sample_id == "sample01"]

out_m <- calculateMoransI(mat1, listw = colGraph(sfe, "visium", sample_id = "sample01"))
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
})

test_that("Correct structure of colGeometryMoransI output", {
    out <- colGeometryMoransI(sfe,
        colGeometryName = "spotPoly",
        colGraphName = "visium", features = "foo",
        sample_id = "sample01"
    )
    fd <- attr(colGeometry(out, "spotPoly", sample_id = "all"), "featureData")
    expect_s4_class(fd, "DataFrame")
    expect_equal(names(fd), c("moran_sample01", "K_sample01"))
    expect_equal(rownames(fd), c("geometry", "foo"))
    expect_true(is.na(fd["geometry", "moran_sample01"]))
    expect_false(is.na(fd["foo", "moran_sample01"]))
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
    expect_equal(rd, out_m)
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
    fd <- attr(colGeometry(out, "spotPoly", sample_id = "all"), "featureData")
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

test_that("DataFrame results for sp.correlogram", {
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
    "localG", "Gi", "E.Gi", "Var.Gi", "Pr(z != E(Gi))",
    "Pr(z != E(Gi)) Sim", "Pr(folded) Sim", "Skewness",
    "Kurtosis", "-log10p Sim", "-log10p_adj Sim"
)
test_that("DataFrame output for localG_perm", {
    out <- calculateUnivariate(mat1,
        listw = colGraph(sfe, "visium", sample_id = "sample01"),
        type = "localG_perm"
    )
    expect_s4_class(out, "DFrame")
    expect_true(all(vapply(out, is.matrix, FUN.VALUE = logical(1))))
    expect_equal(colnames(out[[1]]), names_expect_lg)
    expect_equal(names(out), rownames(mat1))
    expect_equal(nrow(out), ncol(mat1))
})

test_that("DataFrame output for localG, not perm", {
    out <- calculateUnivariate(mat1,
                               listw = colGraph(sfe, "visium", "sample01"),
                               type = "localG")
    expect_s4_class(out, "DFrame")
    expect_true(all(vapply(out, function(o) is.atomic(o) & is.vector(o) &
                               is.numeric(o), FUN.VALUE = logical(1))))
    expect_equal(names(out), rownames(mat1))
    expect_equal(nrow(out), ncol(mat1))
})
