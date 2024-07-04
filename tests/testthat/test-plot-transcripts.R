library(SFEData)
library(SingleCellExperiment)
library(SpatialFeatureExperiment)

expect_ggplot <- function(description, g) {
    expect_s3_class(g, "ggplot")
    expect_error(ggplot_build(g), NA)
}

fp <- tempfile()
fn <- XeniumOutput("v2", file_path = fp)
try(sfe <- readXenium(fn, add_molecules = TRUE))
sfe <- readXenium(fn, add_molecules = TRUE)
bbox <- c(xmin = 0, xmax = 200, ymin = -600, ymax = -400)
test_that("plotGeometry plotting multiple colGeometries", {
    expect_ggplot("Both cells and nuclei",
                  plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"), bbox = bbox))
    expect_ggplot("Both cells and nuclei, dark theme",
                  plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"), bbox = bbox, dark = TRUE))
    expect_ggplot("Both cells and nuclei, with image",
                  plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"), bbox = bbox,
                               dark = TRUE, image_id = "morphology_focus", channel = 3:1))
})

sfe_muscle <- McKellarMuscleData("small")

test_that("plotGeometry plot both colGeometry and annotGeometry", {
    expect_ggplot("Both spotPoly and myofibers",
                  plotGeometry(sfe_muscle, colGeometryName = "spotPoly",
                               annotGeometryName = "myofiber_simplified"))
})

set.seed(29)
genes <- rownames(sfe)[rowData(sfe)$Type == "Gene Expression"]
total <- rowSums(counts(sfe))
probs <- total/sum(total)
genes_use <- sample(rownames(sfe)[rowData(sfe)$Type == "Gene Expression"], 3)

test_that("plotGeometry plotting transcript spots", {
    expect_ggplot("All transcripts regardless of genes",
                  plotGeometry(sfe, rowGeometryName = "txSpots", bbox = bbox, tx_alpha = 0.3))
    expect_ggplot("All transcripts regardless of genes, dark mode",
                  plotGeometry(sfe, rowGeometryName = "txSpots", bbox = bbox,
                               dark = TRUE, tx_alpha = 0.3))
    expect_ggplot("Transcripts of select genes",
                  plotGeometry(sfe, rowGeometryName = "txSpots", gene = genes_use, bbox = bbox,
                               dark = TRUE, tx_alpha = 0.8))
    expect_ggplot("Plot both transcripts and colGeometries",
                  plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"),
                               rowGeometryName = "txSpots", gene = genes_use, bbox = bbox,
                               dark = TRUE, tx_alpha = 0.8))
    expect_ggplot("Plot those geometries with image",
                  plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"),
                               rowGeometryName = "txSpots", gene = genes_use, bbox = bbox,
                               dark = TRUE, tx_alpha = 0.8, image_id = "morphology_focus",
                               channel = 3:1))
})

test_that("Selectively load genes for plotting", {
    skip_if_not(gdalParquetAvailable())
    sfe2 <- readXenium(fn, add_molecules = FALSE, row.names = "symbol")
    expect_ggplot("When no rowGeometry is present",
                  plotGeometry(sfe2, rowGeometryName = "txSpots", gene = genes_use,
                               tx_file = file.path(fn, "tx_spots.parquet")))
    tx <- readSelectTx(file.path(fn, "tx_spots.parquet"), gene_select = head(rownames(sfe2)))
    txSpots(sfe2) <- tx
    expect_ggplot("Load some but not all genes",
                  plotGeometry(sfe2, rowGeometryName = "txSpots", gene = genes_use,
                               tx_file = file.path(fn, "tx_spots.parquet")))
})

test_that("Error when Parquet driver is unavailable", {
    skip_if(gdalParquetAvailable())
    sfe2 <- readXenium(fn, add_molecules = FALSE, row.names = "symbol")
    expect_error(plotGeometry(sfe2, rowGeometryName = "txSpots", gene = genes_use,
                              tx_file = file.path(fn, "tx_spots.parquet")),
                 "GDAL Parquet driver is required")
})

test_that("plotTxBin2D", {
    expect_ggplot("From txSpots of SFE object, square grid",
                  plotTxBin2D(sfe, binwidth = 20))
    expect_ggplot("From txSpots of SFE object, hex grid",
                  plotTxBin2D(sfe, binwidth = 20, hex = TRUE))
    expect_ggplot("A subset of genes, from SFE object",
                  plotTxBin2D(sfe, binwidth = 30, hex = TRUE, gene = genes_use))
    expect_ggplot("From data_dir",
                  plotTxBin2D(data_dir = fn, tech = "Xenium", binwidth = 20, flip = TRUE))
    expect_ggplot("Directly form file",
                  plotTxBin2D(file = file.path(fn, "transcripts.parquet"),
                              binwidth = 20, flip = TRUE,
                              spatialCoordsNames = c("x_location", "y_location"),
                              gene_col = "feature_name"))
    expect_ggplot("For a subset of genes",
                  plotTxBin2D(data_dir = fn, tech = "Xenium", binwidth = 30,
                              flip = TRUE, gene = genes_use))
    # Error when none of sfe, data_dir, or file is specified
    expect_error(plotTxBin2D(), "One of sfe, data_dir, and file must be specified.")
    expect_error(plotTxBin2D(data_dir = fn, tech = "Xenium", binwidth = 30,
                             flip = TRUE, gene = "foobar"),
                 "None of the genes are in the transcript spot file")
    expect_error(plotTxBin2D(sfe, binwidth = 20, gene = "foobar"),
                 "None of the genes are in the SFE object")
})

unlink(fn, recursive = TRUE)
