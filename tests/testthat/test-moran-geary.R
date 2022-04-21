library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(Matrix)
sfe <- readRDS(system.file("testdata/sfe.rds", package = "Voyager"))
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
  colData(sfe)$nCounts <- colSums(mat)
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
  colGeometry(sfe, "spotPoly")$foo <- rnorm(ncol(sfe))
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
