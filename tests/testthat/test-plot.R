# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(spatialLIBD)
library(stringr)
library(vdiffr)
library(scRNAseq)
library(scater)
# Toy example
sfe <- readRDS(system.file("testdata/sfe.rds", package = "Voyager"))
set.seed(29)
colGeometry(sfe, "spotPoly")$foo <- rnorm(ncol(sfe))
mat <- assay(sfe, "counts")
colData(sfe)$nCounts <- colSums(mat)
sfe <- runMoranPlot(sfe, "visium1", c("B", "H"), sample_id = "sample01",
                    exprs_values = "counts")
sfe <- colDataMoranPlot(sfe, "visium1", features = "nCounts",
                        sample_id = "sample01")

# Real dataset, from spatialLIBD
ehub <- ExperimentHub::ExperimentHub()
sce <- fetch_data(type = "sce_example", eh = ehub)
# There's duplicated barcode. Not my fault. Just remove it. It's a toy example.
sce <- sce[,!duplicated(colnames(sce))]
sce <- sce[rowSums(assay(sce, "counts")) > 0,]
spe <- sce_to_spe(sce)
sfe_libd <- toSpatialFeatureExperiment(spe)
sample_use <- "151507"
colGraph(sfe_libd, "visium", sample_id = sample_use) <-
  findVisiumGraph(sfe_libd, sample_use)
# Actually, the toy example only has one gene.
feature_use <- sample(rownames(sfe_libd), 1)
sfe_libd <- runMoranPlot(sfe_libd, "visium", sample_id = sample_use,
                         features = feature_use)
sfe_libd <- colDataMoranPlot(sfe_libd, "visium", features = "sum_umi",
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
                      plotColGraph(sfe, "visium1", colGeometryName = "spotPoly",
                                   sample_id = "sample01"))
})

# Just toy example to test plotting functions
# Need this because the spe toy dataset only has 1 gene
# These functions are not specific to SFE, but I wrote them because I'm not
# satisfied with existing plotting functions.
fluidigm <- ReprocessedFluidigmData(assays = "tophat_counts")
fluidigm <- fluidigm[rowSums(assay(fluidigm, "tophat_counts")) > 0,]
fluidigm <- fluidigm[,colSums(assay(fluidigm, "tophat_counts")) > 0]
tot_counts <- colSums(assay(fluidigm, "tophat_counts"))
logcounts(fluidigm) <- log1p(sweep(assay(fluidigm, "tophat_counts"), 2, tot_counts, "/")*1e5)
fluidigm <- runPCA(fluidigm, ncomponents = 20)

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
