# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(vdiffr)
library(scater)
library(Matrix)
library(ggplot2)
# Toy example
sfe <- readRDS(system.file("extdata/sfe.rds", package = "Voyager"))
sfe <- runUnivariate(sfe,
    type = "moran.plot", c("B", "H"), "visium",
    sample_id = "sample01", exprs_values = "counts"
)

test_that("Everything plotSpatialFeature", {
    expect_doppelganger(
        "Plot gene expression",
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_doppelganger(
        "Plot colData",
        plotSpatialFeature(sfe, "nCounts", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_doppelganger(
        "Plot colGeometry",
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_doppelganger("Plot with annotGeometry", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot"
        )
    })
    expect_doppelganger("Plot with annotGeometry, with new fill scale", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot",
            annot_aes = list(fill = "bar")
        )
    })
    expect_doppelganger("Plot with annotGeometry, colored outlines of polygons", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot",
            annot_aes = list(color = "bar"),
            annot_fixed = list(fill = NA)
        )
    })
    expect_doppelganger("Divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Divergent scale, annot also on divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0, annotGeometryName = "annot",
            annot_aes = list(fill = "bar"),
            annot_divergent = TRUE, annot_diverge_center = 0
        )
    })
    expect_doppelganger("Divergent scale, annot not on divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0, annotGeometryName = "annot",
            annot_aes = list(fill = "bar"),
            annot_divergent = FALSE
        )
    })
    expect_doppelganger("Discrete, represented as point shapes", {
        plotSpatialFeature(sfe, "category", "centroids", "sample01",
            aes_use = "shape", size = 2
        )
    })
    expect_doppelganger("Discrete, represented as color", {
        plotSpatialFeature(sfe, "category", "centroids", "sample01", size = 2)
    })
})

# Real dataset
library(SFEData)
sfe_muscle <- McKellarMuscleData("small")
colGraph(sfe_muscle, "visium") <- findVisiumGraph(sfe_muscle)
sfe_muscle <- logNormCounts(sfe_muscle)
sfe_muscle <- runUnivariate(sfe_muscle, type = "localmoran",
                            c("Myh1", "Myh2"), "visium")

annotGeometry(sfe_muscle, "myofiber_simplified") <-
    sf::st_buffer(annotGeometry(sfe_muscle, "myofiber_simplified"), 0)
annotGraph(sfe_muscle, "poly2nb_myo") <-
    findSpatialNeighbors(sfe_muscle,
        type = "myofiber_simplified", MARGIN = 3,
        method = "poly2nb", zero.policy = TRUE
    )
sfe_muscle <- annotGeometryUnivariate(sfe_muscle, "localmoran",
    features = "area",
    annotGraphName = "poly2nb_myo",
    annotGeometryName = "myofiber_simplified",
    zero.policy = TRUE
)
sfe_muscle <- annotGeometryUnivariate(sfe_muscle, "localG",
                                      features = "area",
                                      annotGraphName = "poly2nb_myo",
                                      annotGeometryName = "myofiber_simplified",
                                      zero.policy = TRUE, include_self = TRUE
)

test_that("Everything plotLocalResult", {
    expect_doppelganger("Plot localmoran Ii for gene", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Plot Ii for gene on top of an annotation", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly",
            annotGeometryName = "tissueBoundary", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Plot another column", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "Z.Ii",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Plot a categorical attribute", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "mean",
            colGeometryName = "spotPoly"
        )
    })
    expect_doppelganger("Plot 2 features", {
        plotLocalResult(sfe_muscle, "localmoran", c("Myh1", "Myh2"),
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Plot Ii for annotGeometry alone", {
        plotLocalResult(sfe_muscle, "localmoran", "area", "Ii",
            annotGeometryName = "myofiber_simplified",
            size = 0.3, color = "cyan", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_doppelganger("Plot a type in annotGeometry but not assay or colData", {
        plotLocalResult(sfe_muscle, "localG", "area",
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
})

feature_use <- "Myh1"
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "moran.plot",
    colGraphName = "visium", features = feature_use,
    exprs_values = "counts"
)
sfe_muscle <- colDataUnivariate(sfe_muscle,
    type = "moran.plot",
    colGraphName = "visium", features = "nCounts"
)
set.seed(29)
colData(sfe_muscle)$GraphBased <- factor(sample(1:5, ncol(sfe_muscle),
    replace = TRUE
),
levels = as.character(1:5)
)

test_that("moranPlot, not filled, no color_by", {
    skip_on_ci()
    expect_warning(
        moranPlot(sfe, "B", "visium1", "sample01"),
        "Too few points"
    )
    expect_doppelganger(
        "moranPlot, not filled",
        moranPlot(sfe_muscle, feature_use, "visium")
    )
    expect_doppelganger(
        "moranPlot, not filled, colData",
        moranPlot(sfe_muscle, "nCounts", "visium")
    )
})

test_that("moranPlot, not filled, with color_by", {
    skip_on_ci()
    expect_doppelganger(
        "moranPlot, not filled, with color_by",
        moranPlot(sfe_muscle, feature_use, "visium",
            color_by = "GraphBased", contour_color = "blue"
        )
    )
})

test_that("moranPlot, filled, no color_by", {
    skip_on_ci()
    expect_doppelganger(
        "moranPlot, filled",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE
        )
    )
})

test_that("moranPlot, filled, with color_by", {
    skip_on_ci()
    expect_doppelganger(
        "moranPlot, filled, with color_by",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE, color_by = "GraphBased"
        )
    )
})

test_that("plotColGraph", {
    expect_doppelganger(
        "plotColGraph toy example",
        plotColGraph(sfe,
            colGraphName = "visium",
            colGeometryName = "spotPoly",
            sample_id = "sample01"
        )
    )
})

sample_use <- "Vis5A"
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use,
    order = 5, zero.policy = TRUE, exprs_values = "counts"
)
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use, order = 5,
    zero.policy = TRUE, method = "corr",
    exprs_values = "counts"
)
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use, order = 5,
    zero.policy = TRUE, method = "C",
    exprs_values = "counts"
)
sfe_muscle <- colDataUnivariate(sfe_muscle,
    type = "sp.correlogram",
    colGraphName = "visium", features = "nCounts",
    sample_id = sample_use, order = 5,
    zero.policy = TRUE
)
test_that("plotCorrelogram", {
    expect_doppelganger(
        "plotCorrelogram, one gene, I",
        plotCorrelogram(sfe_muscle, feature_use, sample_use)
    )
    expect_doppelganger(
        "plotCorrelogram, one gene, corr",
        plotCorrelogram(sfe_muscle, feature_use, sample_use, method = "corr")
    )
    expect_doppelganger(
        "plotCorrelogram, one gene, C",
        plotCorrelogram(sfe_muscle, feature_use, sample_use, method = "C")
    )
    expect_doppelganger(
        "plotCorrelogram, colData, I",
        plotCorrelogram(sfe_muscle, "nCounts", sample_use)
    )
    expect_doppelganger(
        "plotCorrelogram, specify gene and colData, I",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use)
    )
    expect_doppelganger(
        "plotCorrelogram, categorical color_by",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use,
            color_by = c("foo", "bar")
        )
    )
    expect_doppelganger(
        "plotCorrelogram, continuous color_by",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use,
            color_by = 1:2
        ) +
            ggplot2::theme_dark()
    )
})

# Just toy example to test plotting functions
# Need this because the spe toy dataset only has 1 gene
# These functions are not specific to SFE, but I wrote them because I'm not
# satisfied with existing plotting functions.
sfe_muscle <- runPCA(sfe_muscle, ncomponents = 20, BSPARAM = BiocSingular::ExactParam())

test_that("ElbowPlot for PCA", {
    expect_doppelganger("ElbowPlot, default", ElbowPlot(sfe_muscle))
    expect_doppelganger(
        "ElbowPlot, with 10 of the 20 PCs",
        ElbowPlot(sfe_muscle, ndims = 10)
    )
    expect_doppelganger(
        "ElbowPlot, more PCs than available",
        ElbowPlot(sfe_muscle, ndims = 30)
    )
})

test_that("plotDimLoadings for PCA", {
    expect_doppelganger(
        "plotDimLoadings, balanced",
        plotDimLoadings(sfe_muscle, dims = 1:2)
    )
    expect_doppelganger(
        "plotDimLoadings, not balanced",
        plotDimLoadings(sfe_muscle, 1:2, balanced = FALSE)
    )
})

test_that("Everything spatialReducedDim", {
    expect_doppelganger("Plot PCs in space", {
        spatialReducedDim(sfe_muscle, "PCA", 2, "spotPoly",
            annotGeometryName = "tissueBoundary",
            divergent = TRUE, diverge_center = 0
        )
    })
})

test_that("When a gene symbol rowname is not a valid R object name", {
    rownames(sfe_muscle)[1] <- "HLA-foo" # Just toy example
    expect_doppelganger("plotSpatialFeature with illegal gene name",
                        plotSpatialFeature(sfe_muscle, "HLA-foo", "spotPoly",
                                           show_symbol = FALSE))
    sfe_muscle <- runUnivariate(sfe_muscle, "localmoran", "HLA-foo",
                                colGraphName = "visium")
    expect_doppelganger("plotLocalResult with illegal gene name",
                        plotLocalResult(sfe_muscle, "localmoran", "HLA-foo",
                                        colGeometryName = "spotPoly",
                                        show_symbol = FALSE))
})

# scattermore
library(sf)
sfe_cosmx <- HeNSCLCData()
bbox_use <- st_bbox(colGeometry(sfe_cosmx, "centroids"))
annot <- data.frame(x = c(bbox_use[c("xmin", "xmax")]),
                    y = c(bbox_use[c("ymax", "ymin")]),
                    ID = 1)
annot <- df2sf(annot, geometryType = "LINESTRING")
annotGeometry(sfe_cosmx, "foo") <- annot
sfe_cosmx <- sfe_cosmx[, sfe_cosmx$nCounts > 10]
sfe_cosmx <- logNormCounts(sfe_cosmx)
test_that("scattermore in plotSpatialFeature", {
    expect_doppelganger("Plot colData with scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_doppelganger("Plot multiple colData columns", {
        plotSpatialFeature(sfe_cosmx, c("nCounts", "nGenes"),
                           colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_doppelganger("Divergent scale with scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           divergent = TRUE, scattermore = TRUE, size = 0)
    })
    expect_doppelganger("Gene expression", {
        plotSpatialFeature(sfe_cosmx, "KRT19", colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_doppelganger("Also plot annotGeometry", {
        plotSpatialFeature(sfe_cosmx, "KRT19", colGeometryName = "centroids",
                           annotGeometryName = "foo", scattermore = TRUE,
                           size = 0)
    })
    expect_warning(plotSpatialFeature(sfe_cosmx, "nCounts",
                                      colGeometryName = "cellSeg",
                                      scattermore = TRUE),
                   "Using centroids.")
})

localResult(sfe_cosmx, "localG", "KRT19") <- seq_len(ncol(sfe_cosmx))
test_that("Use scattermore in plotLocalResult", {
    expect_doppelganger("scattermore plotLocalResult", {
        plotLocalResult(sfe_cosmx, "localG", "KRT19",
                        colGeometryName = "centroids", scattermore = TRUE,
                        size = 0)
    })
})

fake_pca <- as.matrix(t(logcounts(sfe_cosmx)[c("KRT19", "COL1A1"),]))
colnames(fake_pca) <- c("PC1", "PC2")
reducedDim(sfe_cosmx, "PCA") <- fake_pca
test_that("Use scattermore in spatialReducedDim", {
    expect_doppelganger("scattermore spatialReducedDim", {
        spatialReducedDim(sfe_cosmx, "PCA", 2, scattermore = TRUE)
    })
})

rowData(sfe_cosmx)$is_neg <- grepl("^NegPrb", rownames(sfe_cosmx))
test_that("colData and rowData bin2d", {
    expect_doppelganger("colData bin2d", {
        plotColDataBin2D(sfe_cosmx, "nCounts", "nGenes")
    })
    expect_doppelganger("colData bin2d with hexbin", {
        plotColDataBin2D(sfe_cosmx, "nCounts", "nGenes", hex = TRUE)
    })
    expect_doppelganger("rowData bin2d", {
        plotRowDataBin2D(sfe_cosmx, "means", "vars", bins = 50) +
            scale_x_log10() + scale_y_log10()
    })
    expect_doppelganger("rowData bin2d with subset", {
        plotRowDataBin2D(sfe_cosmx, "means", "vars", subset = "is_neg",
                         name_true = "Counts (negative controls)",
                         name_false = "Counts (real genes)", bins = 50) +
            scale_x_log10() + scale_y_log10()
    })
    expect_doppelganger("rowData bin2d with subset and default legend", {
        plotRowDataBin2D(sfe_cosmx, "means", "vars", subset = "is_neg",
                         bins = 50) +
            scale_x_log10() + scale_y_log10()
    })
})

test_that("colData and rowData histograms", {
    expect_doppelganger("colData histogram, one variable", {
        plotColDataHistogram(sfe_cosmx, "nCounts")
    })
    expect_doppelganger("colData histogram, multiple variables", {
        plotColDataHistogram(sfe_cosmx, c("nCounts", "nGenes"))
    })
    expect_doppelganger("One variable, fill_by", {
        plotRowDataHistogram(sfe_cosmx, "means", fill_by = "is_neg")
    })
    expect_doppelganger("Multiple variables, fill_by", {
        plotRowDataHistogram(sfe_cosmx, c("means", "vars"), fill_by = "is_neg")
    })
    expect_doppelganger("with subset", {
        plotRowDataHistogram(sfe_cosmx, "means", subset = "is_neg")
    })
})

test_that("plotCellBin2D", {
    expect_doppelganger("Cell density, rectangular", {
        plotCellBin2D(sfe_cosmx, bins = 50)
    })
    expect_doppelganger("Cell density, hex", {
        plotCellBin2D(sfe_cosmx, hex = TRUE, bins = 50)
    })
})

sfe_muscle2 <- McKellarMuscleData()
sfe_muscle2 <- sfe_muscle2[, sfe_muscle2$in_tissue]
colGraph(sfe_muscle2, "visium") <- findVisiumGraph(sfe_muscle2)
sfe_muscle2 <- colDataUnivariate(sfe_muscle2, "moran.plot", "nCounts")
test_that("Moran plot bin2d", {
    expect_doppelganger("Moran plot rectangular bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, bins = 30)
    })
    expect_doppelganger("Moran plot hex bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, hex = TRUE, bins = 30)
    })
})
