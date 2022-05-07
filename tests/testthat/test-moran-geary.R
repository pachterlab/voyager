library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(Matrix)
library(bluster)
sfe <- readRDS(system.file("testdata/sfe.rds", package = "Voyager"))
set.seed(29)
mat <- assay(sfe, "counts")
mat1 <- mat[,colData(sfe)$sample_id == "sample01"]

out_m <- calculateMoransI(mat1, listw = colGraph(sfe, "visium1", sample_id = "sample01"))
test_that("Correct structure of calculateMoransI output (matrix)", {
  expect_s4_class(out_m, "DataFrame")
  expect_equal(names(out_m), c("I", "K"))
  expect_true(is.numeric(out_m$I))
  expect_true(is.numeric(out_m$K))
  expect_equal(rownames(out_m), rownames(mat1))
})

test_that("Correct structure of colDataMoransI output", {
  out <- colDataMoransI(sfe, "visium1", "nCounts", sample_id = "sample01")
  expect_s4_class(out, "SpatialFeatureExperiment")
  fd <- attr(colData(out), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), c("MoransI_sample01", "K_sample01"))
  expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
  expect_true(is.na(fd["barcode","MoransI_sample01"]))
  expect_false(is.na(fd["nCounts", "MoransI_sample01"]))
})

test_that("Correct structure of colGeometryMoransI output", {
  out <- colGeometryMoransI(sfe, colGeometryName = "spotPoly",
                            colGraphName = "visium1", features = "foo",
                            sample_id = "sample01")
  fd <- attr(colGeometry(out, "spotPoly"), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), c("MoransI_sample01", "K_sample01"))
  expect_equal(rownames(fd), c("geometry", "foo"))
  expect_true(is.na(fd["geometry","MoransI_sample01"]))
  expect_false(is.na(fd["foo", "MoransI_sample01"]))
})

test_that("Correct structure of calculateMoransI output (SFE)", {
  out <- calculateMoransI(sfe, "visium1", features = rownames(mat1),
                          sample_id = "sample01", exprs_values = "counts")
  expect_s4_class(out, "DataFrame")
  expect_equal(names(out), c("I", "K"))
  expect_true(is.numeric(out$I))
  expect_true(is.numeric(out$K))
  expect_equal(rownames(out), rownames(mat1))
})

test_that("Properly add Moran's I results (no permutation) to SFE rowData", {
  sfe2 <- runMoransI(sfe, "visium1", rownames(mat1), sample_id = "sample01",
                     exprs_values = "counts")
  rd <- rowData(sfe2)
  expect_equal(names(rd), c("MoransI_sample01", "K_sample01"))
  names(rd) <- c("I", "K")
  expect_equal(rd, out_m)
})

test_that("Correct structure of calculateMoranMC output", {
  out <- calculateMoranMC(mat1, listw = colGraph(sfe, "visium1", sample_id = "sample01"),
                          nsim = 10)
  expect_true(is.list(out))
  expect_equal(names(out), rownames(mat1))
  expect_true(all(vapply(out, is, class2 = "htest", FUN.VALUE = logical(1))))
  expect_true(all(vapply(out, is, class2 = "mc.sim", FUN.VALUE = logical(1))))
  i_score <- vapply(out, function(o) unname(o$statistic), FUN.VALUE = numeric(1))
  expect_equal(unname(i_score), out_m$I)
})

names_expect_mc <- c("statistic", "parameter", "p.value", "alternative",
                  "method", "data.name", "res")
names_expect_mc <- paste("MoranMC", names_expect_mc, "sample01", sep = "_")

test_that("Correct structure of colDataMoranMC output", {
  out <- colDataMoranMC(sfe, "visium1", "nCounts", sample_id = "sample01",
                        nsim = 10)
  fd <- attr(colData(out), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), names_expect_mc)
  expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
  expect_true(is.na(fd["barcode","MoranMC_statistic_sample01"]))
  expect_false(is.na(fd["nCounts", "MoranMC_statistic_sample01"]))
})

test_that("Correct structure of colGeometryMoranMC output", {
  out <- colGeometryMoranMC(sfe, "spotPoly", "visium1", "foo", "sample01", 10)
  fd <- attr(colGeometry(out, "spotPoly"), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), names_expect_mc)
  expect_equal(rownames(fd), c("geometry", "foo"))
  expect_true(is.na(fd["geometry", "MoranMC_statistic_sample01"]))
  expect_false(is.na(fd["foo", "MoranMC_statistic_sample01"]))
})

test_that("MoranMC results properly added to rowData", {
  sfe2 <- runMoranMC(sfe, "visium1", rownames(mat1), sample_id = "sample01",
                     exprs_values = "counts", nsim = 10)
  rd <- rowData(sfe2)
  names_expect <- c("statistic", "parameter", "p.value", "alternative",
                    "method", "data.name", "res")
  names_expect <- paste("MoranMC", names_expect, "sample01", sep = "_")
  expect_equal(names(rd), names_expect)
  expect_equal(rd$MoranMC_statistic_sample01, out_m$I)
})

test_that("calculateCorrelogram gives appropriate results (matrix)", {
  out <- calculateCorrelogram(mat1, listw = colGraph(sfe, "visium1", sample_id = "sample01"),
                              order = 2)
  expect_true(all(vapply(out, class, FUN.VALUE = character(1)) == "spcor"))
  expect_equal(names(out), rownames(mat1))
  i1 <- vapply(out, function(o) o$res[1,1], FUN.VALUE = numeric(1))
  expect_equal(unname(i1), out_m$I)
})

test_that("Correct structure of colDataCorrelogram output", {
  out <- colDataCorrelogram(sfe, "visium1", "nCounts", "sample01", order = 2)
  fd <- attr(colData(out), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), "Correlogram_I_sample01")
  expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
  expect_true(is.na(fd["barcode","Correlogram_I_sample01"][[1]]))
  res <- fd$Correlogram_I_sample01[[3]]
  expect_equal(colnames(res), c("I", "expectation", "variance"))
  expect_equal(nrow(res), 2L)
})

test_that("Correct structure of colGeometryCorrelogram output", {
  out <- colGeometryCorrelogram(sfe, "spotPoly", "visium1", "foo", "sample01",
                                order = 2)
  fd <- attr(colGeometry(out, "spotPoly"), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), "Correlogram_I_sample01")
  expect_equal(rownames(fd), c("geometry", "foo"))
  expect_true(is.na(fd["geometry","Correlogram_I_sample01"][[1]]))
  res <- fd$Correlogram_I_sample01[[2]]
  expect_equal(colnames(res), c("I", "expectation", "variance"))
  expect_equal(nrow(res), 2L)
})

test_that("Correlogram results correctly added to rowData", {
  sfe2 <- runCorrelogram(sfe, "visium1", rownames(mat1), sample_id = "sample01",
                         exprs_values = "counts", order = 2)
  rd <- rowData(sfe2)
  expect_equal(names(rd), "Correlogram_I_sample01")
  rdc <- rd$Correlogram_I_sample01
  expect_true(all(vapply(rdc, is.matrix, FUN.VALUE = logical(1))))
  expect_true(all(vapply(rdc, function(r) identical(colnames(r), c("I", "expectation", "variance")),
                         FUN.VALUE = logical(1))))
  i1 <- vapply(rdc, function(r) r[1,1], FUN.VALUE = numeric(1))
  expect_equal(i1, out_m$I)
})

test_that("Correct structure of calculateMoranPlot (matrix)", {
  out <- calculateMoranPlot(mat1, listw = colGraph(sfe, "visium1", sample_id = "sample01"))
  expect_true(all(vapply(out, is.data.frame, FUN.VALUE = logical(1))))
})

names_expect_mp <- c("x", "wx", "is_inf", "labels", "dfb.1_", "dfb.x",
                     "dffit", "cov.r", "cook.d", "hat")

test_that("Correct structure of colDataMoranPlot output", {
  out <- colDataMoranPlot(sfe, "visium1", "nCounts", "sample01")
  fd <- attr(colData(out), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), "MoranPlot_sample01")
  expect_equal(rownames(fd), c("barcode", "sample_id", "nCounts"))
  expect_true(is.na(fd["barcode","MoranPlot_sample01"][[1]]))
  res <- fd$MoranPlot_sample01[[3]]
  expect_equal(colnames(res), names_expect_mp)
  expect_equal(nrow(res), ncol(mat1))
})

test_that("Correct structure of colGeometryMoranPlot output", {
  out <- colGeometryMoranPlot(sfe, "spotPoly", "visium1", "foo", "sample01")
  fd <- attr(colGeometry(out, "spotPoly"), "featureData")
  expect_s4_class(fd, "DataFrame")
  expect_equal(names(fd), "MoranPlot_sample01")
  expect_equal(rownames(fd), c("geometry", "foo"))
  expect_true(is.na(fd["geometry","MoranPlot_sample01"][[1]]))
  res <- fd$MoranPlot_sample01[[2]]
  expect_equal(colnames(res), names_expect_mp)
  expect_equal(nrow(res), ncol(mat1))
})

test_that("Correctly add runMoranPlot output to rowData", {
  sfe2 <- runMoranPlot(sfe, "visium1", rownames(mat1), sample_id = "sample01",
                       exprs_values = "counts")
  rd <- rowData(sfe2)
  expect_equal(names(rd), "MoranPlot_sample01")
  rdc <- rd$MoranPlot_sample01
  expect_true(all(vapply(rdc, is.data.frame, FUN.VALUE = logical(1))))
  expect_true(all(vapply(rdc, function(r) identical(names(r), names_expect_mp),
                         FUN.VALUE = logical(1))))
  expect_true(all(vapply(rdc, nrow, FUN.VALUE = numeric(1)) == ncol(mat1)))
})

# Should have passed the above unit tests for this to work
sfe <- runMoranPlot(sfe, "visium1", c("B", "H"), sample_id = "sample01",
                    exprs_values = "counts")
sfe <- colDataMoranPlot(sfe, "visium1", "nCounts", "sample01")
sfe <- colGeometryMoranPlot(sfe, "spotPoly", "visium1", features = "foo",
                            sample_id = "sample01")

test_that("Moran plot clustering gives right results for gene expression", {
  out <- clusterMoranPlot(sfe, c("B", "H"), KmeansParam(2),
                          sample_id = "sample01")
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("B", "H"))
  expect_true(all(vapply(out, is.factor, FUN.VALUE = logical(1))))
  expect_equal(nrow(out), sum(colData(sfe)$sample_id == "sample01"))
  expect_equal(rownames(out), colnames(sfe)[colData(sfe)$sample_id == "sample01"])
})

test_that("Warning when some of the requested features don't have Moran plot", {
  expect_warning(out <- clusterMoranPlot(sfe, c("B", "H", "L"), KmeansParam(2),
                                         sample_id = "sample01"),
                 "don't have the requested metadata")
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("B", "H"))
})

test_that("Error when none of the features have Moran plot", {
  expect_error(clusterMoranPlot(sfe, c("Q", "L"), KmeansParam(2), "sample01"),
               "None of the features")
})

test_that("Correct results when doing both gene expression and colData", {
  out <- clusterMoranPlot(sfe, c("nCounts", "B", "H"), KmeansParam(2),
                          sample_id = "sample01")
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("B", "H", "nCounts"))
  expect_true(all(vapply(out, is.factor, FUN.VALUE = logical(1))))
  expect_equal(nrow(out), sum(colData(sfe)$sample_id == "sample01"))
  expect_equal(rownames(out), colnames(sfe)[colData(sfe)$sample_id == "sample01"])
})

test_that("Correct Moran plot cluster results for colGeometry", {
  out <- clusterMoranPlot(sfe, "foo", KmeansParam(2), sample_id = "sample01",
                          colGeometryName = "spotPoly")
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), "foo")
  expect_s3_class(out$foo, "factor")
})

test_that("Error when the MoranPlot_sample01 column is absent", {
  rowData(sfe)$MoranPlot_sample01 <- NULL
  expect_error(clusterMoranPlot(sfe, c("Q", "L"), KmeansParam(2), "sample01"),
               "None of the features")
})
