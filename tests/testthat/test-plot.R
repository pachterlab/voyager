# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(spatialLIBD)
library(stringr)
library(vdiffr)
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
spe <- sce_to_spe(sce)
# There's duplicated barcode. Not my fault. Just remove it. It's a toy example.
spe <- spe[,!duplicated(colnames(spe))]
spe <- spe[rowSums(assay(spe, "counts")) > 0,]
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
