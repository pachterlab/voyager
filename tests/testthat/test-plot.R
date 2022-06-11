# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(spatialLIBD)
library(vdiffr)
library(scRNAseq)
library(scater)
library(Matrix)
# Toy example
sfe <- readRDS(system.file("testdata/sfe.rds", package = "Voyager"))
sfe <- runMoranPlot(sfe, c("B", "H"), "visium", sample_id = "sample01",
                    exprs_values = "counts")

test_that("Everything plotSpatialFeature", {
  #testthat::skip("Skipping plots that require sf")
  expect_doppelganger("Plot gene expression",
                      plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
                                         exprs_values = "counts"))
  expect_doppelganger("Plot colData",
                      plotSpatialFeature(sfe, "nCounts", "spotPoly", "sample01",
                                         exprs_values = "counts"))
  expect_doppelganger("Plot colGeometry",
                      plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                                         exprs_values = "counts"))
  expect_doppelganger("Plot with annotGeometry", {
    plotSpatialFeature(sfe, "H", "spotPoly", "sample01", exprs_values = "counts",
                       annotGeometryName = "annot")
  })
  expect_doppelganger("Plot with annotGeometry, with new fill scale", {
    plotSpatialFeature(sfe, "H", "spotPoly", "sample01", exprs_values = "counts",
                       annotGeometryName = "annot",
                       annot_aes = list(fill = "bar"))
  })
  expect_doppelganger("Plot with annotGeometry, colored outlines of polygons", {
    plotSpatialFeature(sfe, "H", "spotPoly", "sample01", exprs_values = "counts",
                       annotGeometryName = "annot",
                       annot_aes = list(color = "bar"),
                       annot_fixed = list(fill = NA))
  })
  expect_doppelganger("Divergent scale", {
    plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                       exprs_values = "counts", divergent = TRUE,
                       diverge_center = 0)
  })
  expect_doppelganger("Divergent scale, annot also on divergent scale", {
    plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                       exprs_values = "counts", divergent = TRUE,
                       diverge_center = 0, annotGeometryName = "annot",
                       annot_aes = list(fill = "bar"),
                       annot_divergent = TRUE, annot_diverge_center = 0)
  })
  expect_doppelganger("Divergent scale, annot not on divergent scale", {
    plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                       exprs_values = "counts", divergent = TRUE,
                       diverge_center = 0, annotGeometryName = "annot",
                       annot_aes = list(fill = "bar"),
                       annot_divergent = FALSE)
  })
  expect_doppelganger("Discrete, represented as point shapes", {
    plotSpatialFeature(sfe, "category", "centroids", "sample01",
                       aes_use = "shape", size = 2)
  })
  expect_doppelganger("Discrete, represented as color", {
    plotSpatialFeature(sfe, "category", "centroids", "sample01", size = 2)
  })
})

# Real dataset, from spatialLIBD
ehub <- ExperimentHub::ExperimentHub()
sce <- fetch_data(type = "sce_example", eh = ehub)
# There's duplicated barcode. Not my fault. Just remove it. It's a toy example.
sce <- sce[,!duplicated(colnames(sce))]
sce <- sce[Matrix::rowSums(assay(sce, "counts")) > 0,]
spe <- sce_to_spe(sce)
sfe_libd <- toSpatialFeatureExperiment(spe)
sample_use <- "151507"
colGraph(sfe_libd, "visium", sample_id = sample_use) <-
  findVisiumGraph(sfe_libd, sample_use)
# Actually, the toy example only has one gene.
feature_use <- sample(rownames(sfe_libd), 1)
sfe_libd <- runMoranPlot(sfe_libd, colGraphName = "visium", sample_id = sample_use,
                         features = feature_use)
sfe_libd <- colDataMoranPlot(sfe_libd, colGraphName = "visium", features = "sum_umi",
                             sample_id = sample_use)
colData(sfe_libd)$GraphBased <- factor(colData(sfe_libd)$GraphBased)

test_that("moranPlot, not filled, no color_by", {
  expect_warning(moranPlot(sfe, "B", "visium1", "sample01"),
                 "Too few points")
  expect_doppelganger("moranPlot, not filled",
                      moranPlot(sfe_libd, feature_use, "visium", sample_use))
  expect_doppelganger("moranPlot, not filled, colData",
                      moranPlot(sfe_libd, "sum_umi", "visium", sample_use))
})

test_that("moranPlot, not filled, with color_by", {
  expect_doppelganger("moranPlot, not filled, with color_by",
                      moranPlot(sfe_libd, feature_use, "visium", sample_use,
                                color_by = "GraphBased", contour_color = "blue"))
})

test_that("moranPlot, filled, no color_by", {
  expect_doppelganger("moranPlot, filled",
                      moranPlot(sfe_libd, feature_use, "visium", sample_use,
                                filled = TRUE))
})

test_that("moranPlot, filled, with color_by", {
  expect_doppelganger("moranPlot, filled, with color_by",
                      moranPlot(sfe_libd, feature_use, "visium", sample_use,
                                filled = TRUE, color_by = "GraphBased"))
})

test_that("plotColGraph", {
  expect_doppelganger("plotColGraph toy example",
                      plotColGraph(sfe, colGraphName = "visium1",
                                   colGeometryName = "spotPoly",
                                   sample_id = "sample01"))
})

# Some spots don't have such high order of neighbors
sfe_libd <- runCorrelogram(sfe_libd, features = feature_use,
                           colGraphName = "visium", sample_id = sample_use,
                           order = 5, zero.policy = TRUE)
sfe_libd <- runCorrelogram(sfe_libd, features = feature_use,
                           colGraphName = "visium", sample_id = sample_use,
                           order = 5, zero.policy = TRUE,
                           method = "corr")
sfe_libd <- runCorrelogram(sfe_libd, features = feature_use,
                           colGraphName = "visium", sample_id = sample_use,
                           order = 5, zero.policy = TRUE,
                           method = "C")
sfe_libd <- colDataCorrelogram(sfe_libd, colGraphName = "visium", features = "sum_umi",
                               sample_id = sample_use, order = 5,
                               zero.policy = TRUE)
test_that("plotCorrelogram", {
  expect_doppelganger("plotCorrelogram, one gene, I",
                      plotCorrelogram(sfe_libd, feature_use, sample_use))
  expect_doppelganger("plotCorrelogram, one gene, corr",
                      plotCorrelogram(sfe_libd, feature_use, sample_use, method = "corr"))
  expect_doppelganger("plotCorrelogram, one gene, C",
                      plotCorrelogram(sfe_libd, feature_use, sample_use, method = "C"))
  expect_doppelganger("plotCorrelogram, colData, I",
                      plotCorrelogram(sfe_libd, "sum_umi", sample_use))
  expect_doppelganger("plotCorrelogram, specify gene and colData, I",
                      plotCorrelogram(sfe_libd, c(feature_use, "sum_umi"), sample_use))
  expect_doppelganger("plotCorrelogram, categorical color_by",
                      plotCorrelogram(sfe_libd, c(feature_use, "sum_umi"), sample_use,
                                      color_by = c("foo", "bar")))
  expect_doppelganger("plotCorrelogram, continuous color_by",
                      plotCorrelogram(sfe_libd, c(feature_use, "sum_umi"), sample_use,
                                      color_by = 1:2) +
                        ggplot2::theme_dark())
})

# Just toy example to test plotting functions
# Need this because the spe toy dataset only has 1 gene
# These functions are not specific to SFE, but I wrote them because I'm not
# satisfied with existing plotting functions.
fluidigm <- ReprocessedFluidigmData(assays = "tophat_counts")
fluidigm <- fluidigm[Matrix::rowSums(assay(fluidigm, "tophat_counts")) > 0,]
fluidigm <- fluidigm[,Matrix::colSums(assay(fluidigm, "tophat_counts")) > 0]
tot_counts <- colSums(assay(fluidigm, "tophat_counts"))
logcounts(fluidigm) <- log1p(sweep(assay(fluidigm, "tophat_counts"), 2, tot_counts, "/")*1e5)
fluidigm <- runPCA(fluidigm, ncomponents = 20, BSPARAM = BiocSingular::ExactParam())

test_that("ElbowPlot for PCA", {
  expect_doppelganger("ElbowPlot, default", ElbowPlot(fluidigm))
  expect_doppelganger("ElbowPlot, with 10 of the 20 PCs",
                      ElbowPlot(fluidigm, ndims = 10))
  expect_doppelganger("ElbowPlot, more PCs than available",
                      ElbowPlot(fluidigm, ndims = 30))
})

test_that("plotDimLoadings for PCA", {
  expect_doppelganger("plotDimLoadings, balanced",
                      plotDimLoadings(fluidigm, dims = 1:2))
  expect_doppelganger("plotDimLoadings, not balanced",
                      plotDimLoadings(fluidigm, 1:2, balanced = FALSE))
})
