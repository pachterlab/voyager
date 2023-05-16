# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(vdiffr)
library(scater)
library(Matrix)
library(ggplot2)
library(scran)

expect_ggplot <- function(description, g) {
    expect_s3_class(g, "ggplot")
    expect_error(ggplot_build(g), NA)
}

test_that("Divergent palette beginning and end", {
    expect_equal(getDivergeRange(1:8, diverge_center = 5), c(0, 0.875))
    expect_equal(getDivergeRange(1:8, diverge_center = 4), c(0.125, 1))
    expect_equal(getDivergeRange(1:8, diverge_center = 0), c(9/16, 1))
    expect_equal(getDivergeRange(1:8, diverge_center = 9), c(0, 8/18))
})

# Toy example
sfe <- readRDS(system.file("extdata/sfe.rds", package = "Voyager"))
sfe <- runUnivariate(sfe,
    type = "moran.plot", c("B", "H"), "visium",
    sample_id = "sample01", exprs_values = "counts"
)

test_that("Everything plotSpatialFeature", {
    expect_ggplot(
        "Plot gene expression",
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_ggplot("Plot gene expression, dark theme",
                  plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
                                     exprs_values = "counts", dark = TRUE
                  ))
    expect_ggplot(
        "Plot colData",
        plotSpatialFeature(sfe, "nCounts", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_ggplot(
        "Plot colGeometry",
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts"
        )
    )
    expect_ggplot("Plot with annotGeometry", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot"
        )
    })
    expect_ggplot("Plot with annotGeometry, with new fill scale", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot",
            annot_aes = list(fill = "bar")
        )
    })
    expect_ggplot("Plot with annotGeometry, with new fill scale, dark theme", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
                           exprs_values = "counts",
                           annotGeometryName = "annot",
                           annot_aes = list(fill = "bar")
        )
    })
    expect_ggplot("Plot with annotGeometry, colored outlines of polygons", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
            exprs_values = "counts",
            annotGeometryName = "annot",
            annot_aes = list(color = "bar"),
            annot_fixed = list(fill = NA)
        )
    })
    expect_ggplot("Plot with annotGeometry, colored outlines of polygons, dark theme", {
        plotSpatialFeature(sfe, "H", "spotPoly", "sample01",
                           exprs_values = "counts",
                           annotGeometryName = "annot",
                           annot_aes = list(color = "bar"),
                           annot_fixed = list(fill = NA), dark = TRUE
        )
    })
    expect_ggplot("Divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_ggplot("Divergent scale, dark theme", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                           exprs_values = "counts", divergent = TRUE,
                           diverge_center = 0, dark = TRUE
        )
    })
    expect_ggplot("Divergent scale, annot also on divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0, annotGeometryName = "annot",
            annot_aes = list(fill = "bar"),
            annot_divergent = TRUE, annot_diverge_center = 0
        )
    })
    expect_ggplot("Divergent scale, annot also on divergent scale, dark theme", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
                           exprs_values = "counts", divergent = TRUE,
                           diverge_center = 0, annotGeometryName = "annot",
                           annot_aes = list(fill = "bar"),
                           annot_divergent = TRUE, annot_diverge_center = 0,
                           dark = TRUE
        )
    })
    expect_ggplot("Divergent scale, annot not on divergent scale", {
        plotSpatialFeature(sfe, "foo", "spotPoly", "sample01",
            exprs_values = "counts", divergent = TRUE,
            diverge_center = 0, annotGeometryName = "annot",
            annot_aes = list(fill = "bar"),
            annot_divergent = FALSE
        )
    })
    expect_ggplot("Discrete, represented as point shapes", {
        plotSpatialFeature(sfe, "category", "centroids", "sample01",
            aes_use = "shape", size = 2
        )
    })
    expect_ggplot("Discrete, represented as color", {
        plotSpatialFeature(sfe, "category", "centroids", "sample01", size = 2)
    })
})

# Real dataset
library(SFEData)
sfe_muscle <- McKellarMuscleData("small")
colGraph(sfe_muscle, "visium") <- findVisiumGraph(sfe_muscle)
sfe_muscle <- logNormCounts(sfe_muscle)
sfe_muscle <- runUnivariate(sfe_muscle, type = "localmoran",
                            c("Myh1", "Myh2"), "visium",
                            swap_rownames = "symbol")

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
    expect_ggplot("Plot localmoran Ii for gene", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_ggplot("Plot localmoran Ii for gene, dark theme", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
                        colGeometryName = "spotPoly", divergent = TRUE,
                        diverge_center = 0, swap_rownames = "symbol",
                        dark = TRUE
        )
    })
    expect_ggplot("Plot Ii for gene on top of an annotation", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly",
            annotGeometryName = "tissueBoundary", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_ggplot("Plot another column", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "Z.Ii",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_ggplot("Plot a categorical attribute", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "mean",
            colGeometryName = "spotPoly", swap_rownames = "symbol"
        )
    })
    expect_ggplot("Plot 2 features", {
        plotLocalResult(sfe_muscle, "localmoran", c("Myh1", "Myh2"),
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_ggplot("Plot Ii for annotGeometry alone", {
        plotLocalResult(sfe_muscle, "localmoran", "area", "Ii",
            annotGeometryName = "myofiber_simplified",
            linewidth = 0.3, color = "cyan", divergent = TRUE,
            diverge_center = 0
        )
    })
    expect_ggplot("Plot Ii for annotGeometry alone, dark theme", {
        plotLocalResult(sfe_muscle, "localmoran", "area", "Ii",
                        annotGeometryName = "myofiber_simplified",
                        linewidth = 0.3, color = "navy", divergent = TRUE,
                        diverge_center = 0, dark = TRUE
        )
    })
    expect_ggplot("Plot a type in annotGeometry but not assay or colData", {
        plotLocalResult(sfe_muscle, "localG", "area",
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
})

feature_use <- "Myh1"
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "moran.plot",
    colGraphName = "visium", features = feature_use,
    exprs_values = "counts", swap_rownames = "symbol"
)
sfe_muscle <- colDataUnivariate(sfe_muscle,
    type = "moran.plot",
    colGraphName = "visium", features = "nCounts"
)

set.seed(29)
colData(sfe_muscle)$GraphBased <- factor(sample(1:5, ncol(sfe_muscle),
    replace = TRUE), levels = as.character(1:5))

test_that("moranPlot, not filled, no color_by", {
    expect_warning(
        moranPlot(sfe, "B", "visium", "sample01"),
        "Too few points"
    )
    expect_ggplot("",
        moranPlot(sfe_muscle, feature_use, "visium", swap_rownames = "symbol")
    )
    expect_ggplot("",
        moranPlot(sfe_muscle, "nCounts", "visium")
    )
})

test_that("moranPlot, not filled, with color_by", {
    expect_ggplot("",
        moranPlot(sfe_muscle, feature_use, "visium",
            color_by = "GraphBased", contour_color = "blue",
            swap_rownames = "symbol"
        )
    )
})

test_that("moranPlot, filled, no color_by", {
    expect_ggplot("",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE, swap_rownames = "symbol"
        )
    )
})

test_that("moranPlot, filled, with color_by", {
    expect_ggplot("",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE, color_by = "GraphBased", swap_rownames = "symbol"
        )
    )
})

sfe_muscle <- annotGeometryUnivariate(sfe_muscle, "moran.plot",
                                      features = "area",
                                      annotGraphName = "poly2nb_myo",
                                      annotGeometryName = "myofiber_simplified",
                                      zero.policy = TRUE
)

test_that("moranPlot, with singletons", {
    expect_ggplot("not filled", {
        moranPlot(sfe_muscle, "area", "poly2nb_myo", annotGeometryName = "myofiber_simplified")
    })
    expect_ggplot("filled", {
        moranPlot(sfe_muscle, "area", "poly2nb_myo",
                  annotGeometryName = "myofiber_simplified", filled = TRUE)
    })
})

test_that("plot graphs", {
    expect_ggplot(
        "plotColGraph",
        plotColGraph(sfe_muscle, colGraphName = "visium",
                     colGeometryName = "spotPoly")
    )
    expect_ggplot("plotAnnotGraph",
                        plotAnnotGraph(sfe_muscle, "poly2nb_myo", "myofiber_simplified"))
})

sample_use <- "Vis5A"
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use, order = 5, zero.policy = TRUE,
    exprs_values = "counts", swap_rownames = "symbol"
)
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use, order = 5,
    zero.policy = TRUE, method = "corr",
    exprs_values = "counts", swap_rownames = "symbol"
)
sfe_muscle <- runUnivariate(sfe_muscle,
    type = "sp.correlogram",
    features = feature_use, colGraphName = "visium",
    sample_id = sample_use, order = 5,
    zero.policy = TRUE, method = "C",
    exprs_values = "counts", swap_rownames = "symbol"
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
        plotCorrelogram(sfe_muscle, feature_use, sample_use,
                        swap_rownames = "symbol")
    )
    expect_doppelganger(
        "plotCorrelogram, one gene, corr",
        plotCorrelogram(sfe_muscle, feature_use, sample_use,
                        method = "corr", swap_rownames = "symbol")
    )
    expect_doppelganger(
        "plotCorrelogram, one gene, C",
        plotCorrelogram(sfe_muscle, feature_use, sample_use,
                        method = "C", swap_rownames = "symbol")
    )
    expect_doppelganger(
        "plotCorrelogram, colData, I",
        plotCorrelogram(sfe_muscle, "nCounts", sample_use)
    )
    expect_doppelganger(
        "plotCorrelogram, specify gene and colData, I",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use,
                        swap_rownames = "symbol")
    )
    expect_doppelganger(
        "plotCorrelogram, categorical color_by",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use,
            color_by = c("foo", "bar"), swap_rownames = "symbol"
        )
    )
    expect_doppelganger(
        "plotCorrelogram, continuous color_by",
        plotCorrelogram(sfe_muscle, c(feature_use, "nCounts"), sample_use,
            color_by = 1:2, swap_rownames = "symbol"
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

sfe1 <- McKellarMuscleData("small")
sfe1 <- sfe1[,sfe1$in_tissue]
sfe1 <- logNormCounts(sfe1)
inds <- order(Matrix::rowSums(logcounts(sfe1)), decreasing = TRUE)[1:50]
sfe2 <- McKellarMuscleData("small2")
sfe2 <- sfe2[,sfe2$in_tissue]
sfe2 <- logNormCounts(sfe2)
sfe3 <- SpatialFeatureExperiment::cbind(sfe1, sfe2)
colGraphs(sfe3, sample_id = "all", name = "visium") <-
    findVisiumGraph(sfe3, "all")

sfe3 <- runMultivariate(sfe3, "multispati", colGraphName = "visium",
                        subset_row = inds, sample_action = "separate",
                        nfposi = 10, nfnega = 10)
test_that("Plot elbow plots for multiple samples", {
    expect_doppelganger("Not facetting multispati elbow, positive only", {
        ElbowPlot(sfe3, facet = FALSE, reduction = "multispati")
    })
    expect_doppelganger("Not facetting multispati elbow, negative only", {
        ElbowPlot(sfe3, facet = FALSE, reduction = "multispati",
                  ndims = 0, nfnega = 10)
    })
    expect_doppelganger("Both positive and negative", {
        ElbowPlot(sfe3, facet = FALSE, reduction = "multispati",
                  ndims = 10, nfnega = 10)
    })
    expect_doppelganger("Facetting by sample", {
        ElbowPlot(sfe3, facet = TRUE, reduction = "multispati",
                  ndims = 10, nfnega = 10)
    })
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
    expect_doppelganger("Change the number of columns",
                        plotDimLoadings(sfe_muscle, 1:2, balanced = TRUE,
                                        ncol = 1))
})

test_that("plotDimLoadings for multiple samples", {
    expect_doppelganger("Multispati loadings, 2 samples", {
        plotDimLoadings(sfe3, c(1:2, 19:20), swap_rownames = "symbol",
                        reduction = "multispati")
    })
})

test_that("Everything spatialReducedDim", {
    expect_ggplot("Plot PCs in space", {
        spatialReducedDim(sfe_muscle, "PCA", 2, colGeometryName = "spotPoly",
            annotGeometryName = "tissueBoundary",
            divergent = TRUE, diverge_center = 0
        )
    })
    expect_ggplot("Plot only one component", {
        spatialReducedDim(sfe_muscle, "PCA", components = 2,
                          colGeometryName = "spotPoly",
                          divergent = TRUE, diverge_center = 0)
    })
})

test_that("When a gene symbol rowname is not a valid R object name", {
    rownames(sfe_muscle)[1] <- "HLA-foo" # Just toy example
    expect_ggplot("plotSpatialFeature with illegal gene name",
                        plotSpatialFeature(sfe_muscle, "HLA-foo", "spotPoly"))
    expect_ggplot("plotSpatialFeature with both illegal and legal names", {
        plotSpatialFeature(sfe_muscle, c("HLA-foo", rownames(sfe_muscle)[3]),
                           "spotPoly")
    })
    rownames(sfe_muscle)[2] <- "HLA-bar"
    expect_ggplot("plotSpatialFeature with multiple illegal names", {
        plotSpatialFeature(sfe_muscle, c("HLA-foo", "HLA-bar"),
                           "spotPoly")
    })
    sfe_muscle <- runUnivariate(sfe_muscle, "localmoran", "HLA-foo",
                                colGraphName = "visium")
    expect_ggplot("plotLocalResult with illegal gene name",
                        plotLocalResult(sfe_muscle, "localmoran", "HLA-foo",
                                        colGeometryName = "spotPoly"))
})

inds <- c(1,3,4,5)
sfe3 <- runUnivariate(sfe3, "sp.correlogram", features = rownames(sfe3)[inds],
                      order = 5)
test_that("Plot correlogram for multiple samples", {
    expect_doppelganger("plotCorrelogram, multiple samples", {
        plotCorrelogram(sfe3, rownames(sfe3)[inds], sample_id = "all")
    })
    expect_doppelganger("Cluster by gene instead", {
        plotCorrelogram(sfe3, rownames(sfe3)[inds], sample_id = "all",
                        facet_by = "features")
    })
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
    expect_ggplot("Plot colData with scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_ggplot("Plot multiple colData columns", {
        plotSpatialFeature(sfe_cosmx, c("nCounts", "nGenes"),
                           colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_ggplot("Divergent scale with scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           divergent = TRUE, scattermore = TRUE, size = 0)
    })
    expect_ggplot("Gene expression", {
        plotSpatialFeature(sfe_cosmx, "KRT19", colGeometryName = "centroids",
                           scattermore = TRUE, size = 0)
    })
    expect_ggplot("Also plot annotGeometry", {
        plotSpatialFeature(sfe_cosmx, "KRT19", colGeometryName = "centroids",
                           annotGeometryName = "foo", scattermore = TRUE,
                           size = 0)
    })
    expect_message(plotSpatialFeature(sfe_cosmx, "nCounts",
                                      colGeometryName = "cellSeg",
                                      scattermore = TRUE),
                   "Using centroids.")
})

localResult(sfe_cosmx, "localG", "KRT19") <- seq_len(ncol(sfe_cosmx))
test_that("Use scattermore in plotLocalResult", {
    expect_ggplot("scattermore plotLocalResult", {
        plotLocalResult(sfe_cosmx, "localG", "KRT19",
                        colGeometryName = "centroids", scattermore = TRUE,
                        size = 0)
    })
})

fake_pca <- as.matrix(t(logcounts(sfe_cosmx)[c("KRT19", "COL1A1"),]))
colnames(fake_pca) <- c("PC1", "PC2")
reducedDim(sfe_cosmx, "PCA") <- fake_pca
test_that("Use scattermore in spatialReducedDim", {
    expect_ggplot("scattermore spatialReducedDim", {
        spatialReducedDim(sfe_cosmx, "PCA", 2, scattermore = TRUE)
    })
})

rowData(sfe_cosmx)$is_neg <- grepl("^NegPrb", rownames(sfe_cosmx))

inds <- spatialCoords(sfe_cosmx)[,1] > 20000
sfe_cosmx2a <- sfe_cosmx[,!inds]
sfe_cosmx2b <- sfe_cosmx[,inds]
colData(sfe_cosmx2b)$sample_id <- "sample02"
names(int_metadata(sfe_cosmx2b)$spatialGraphs) <- "sample02"
sfe_cosmx2 <- SpatialFeatureExperiment::cbind(sfe_cosmx2a, sfe_cosmx2b)
sfe_cosmx2 <- removeEmptySpace(sfe_cosmx2)

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
    expect_doppelganger("colData bin2d for multiple samples", {
        plotColDataBin2D(sfe_cosmx2, "nCounts", "nGenes",
                         facet_by = "sample_id")
    })
    expect_warning(plotColDataBin2D(sfe_cosmx, "nCounts", "nGenes",
                                    facet_by = "nCounts"),
                   "Not facetting")
    expect_warning(plotColDataBin2D(sce = sfe_cosmx, x = "nCounts",
                                    y = "nGenes"),
                   "deprecated")
})

test_that("colData and rowData histograms", {
    expect_doppelganger("colData histogram, one variable", {
        plotColDataHistogram(sfe_cosmx, "nCounts")
    })
    expect_doppelganger("colData histogram, multiple variables", {
        plotColDataHistogram(sfe_cosmx, c("nCounts", "nGenes"))
    })
    expect_doppelganger("One variable, fill_by", {
        plotRowDataHistogram(sfe_cosmx, "means", fill_by = "is_neg",
                             position = "stack")
    })
    expect_doppelganger("Multiple variables, fill_by", {
        plotRowDataHistogram(sfe_cosmx, c("means", "vars"), fill_by = "is_neg",
                             position = "stack")
    })
    expect_doppelganger("with subset", {
        plotRowDataHistogram(sfe_cosmx, "means", subset = "is_neg")
    })
    expect_doppelganger("One variable, facetting", {
        plotColDataHistogram(sfe_cosmx2, "nCounts", facet_by = "sample_id")
    })
    expect_doppelganger("Multiple variables, facetting", {
        plotColDataHistogram(sfe_cosmx2, c("nCounts", "nGenes"),
                             facet_by = "sample_id")
    })
    expect_doppelganger("Illegal symbol", {
        sfe_cosmx$`n-counts` <- sfe_cosmx$nCounts
        plotColDataHistogram(sfe_cosmx, "n-counts")
    })
    expect_warning(plotColDataHistogram(sfe = sfe_cosmx, feature = "nCounts"),
                   "deprecated")
})

test_that("colData and rowData freqpoly", {
    expect_doppelganger("colData freqpoly, one variable", {
        plotColDataFreqpoly(sfe_cosmx, "nCounts")
    })
    expect_doppelganger("colData freqpoly, multiple variables", {
        plotColDataFreqpoly(sfe_cosmx, c("nCounts", "nGenes"))
    })
    expect_doppelganger("One variable, color_by", {
        plotRowDataFreqpoly(sfe_cosmx, "means", color_by = "is_neg")
    })
    expect_doppelganger("Multiple variables, color_by", {
        plotRowDataFreqpoly(sfe_cosmx, c("means", "vars"), color_by = "is_neg")
    })
    expect_doppelganger("with subset, freqpoly", {
        plotRowDataFreqpoly(sfe_cosmx, "means", subset = "is_neg")
    })
    expect_doppelganger("Illegal symbol freqpoly", {
        sfe_cosmx$`n-counts` <- sfe_cosmx$nCounts
        plotColDataFreqpoly(sfe_cosmx, "n-counts")
    })
})

test_that("plotCellBin2D", {
    expect_doppelganger("Cell density, rectangular", {
        plotCellBin2D(sfe_cosmx, bins = 50)
    })
    expect_doppelganger("Cell density, hex", {
        plotCellBin2D(sfe_cosmx, hex = TRUE, bins = 50)
    })
    expect_doppelganger("Multiple samples", {
        plotCellBin2D(sfe_cosmx2, bins = 50)
    })
})

test_that("Binning values", {
    expect_ggplot("Bin and summarize feature", {
        plotSpatialFeature(sfe_cosmx, features = "nCounts",
                           colGeometryName = "centroids", bins = 50)
    })
    expect_ggplot("Bin and summarize dimension reduction values", {
        spatialReducedDim(sfe_cosmx, dimred = "PCA", ncomponents = 1,
                          colGeometryName = "centroids", bins = 50)
    })
    expect_ggplot("Bin and summarize local results", {
        plotLocalResult(sfe_cosmx, "localG", "KRT19",
                        colGeometryName = "centroids", bins = 50)
    })

})
# Using bbox
sfe1 <- McKellarMuscleData("small")
sfe2 <- McKellarMuscleData("small2")
sfe <- SpatialFeatureExperiment::cbind(sfe1, sfe2)
sfe <- removeEmptySpace(sfe)

annotGeometry(sfe, "myofiber_simplified", sample_id = "all", translate = FALSE) <-
    sf::st_buffer(annotGeometry(sfe, "myofiber_simplified", sample_id = "all"), 0)

bbox <- c(xmin = 5000, ymin = 175000, xmax = 6000, ymax = 176000)
bbox_large <- c(xmin = 10000, xmax = 20000, ymin = 160000, ymax = 170000)
bbox2 <- c(xmin = 5500, ymin = 13500, xmax = 6500, ymax = 14500)
bbox_2s1 <- c(xmin = 500, xmax = 1500, ymin = 500, ymax = 1500)
bbox_2s2 <- c(xmin = 0, xmax = 1000, ymin = 1000, ymax = 2000)
bbox_2s <- cbind(bbox_2s1, bbox_2s2)
colnames(bbox_2s) <- c("Vis5A", "sample02")
test_that("Using bbox with plotSpatialFeature", {
    cat("beginning bbox tests")
    # One sample
    expect_ggplot("Only plotting colGeometry", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "cellSeg",
                           bbox = bbox)
    })
    expect_ggplot("With scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           bbox = bbox_large, scattermore = TRUE, pointsize = 1)
    })
    expect_ggplot("Both colGeometry and annotGeometry", {
        plotSpatialFeature(sfe_muscle, "nCounts",
                           annotGeometryName = "myofiber_simplified",
                           annot_aes = list(fill = "area"), bbox = bbox2)
    })
    # Two samples
    expect_ggplot("Two samples, only plotting colGeometry, same bbox", {
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s1)
    })
    expect_ggplot("Two samples, with annotGeometry, same bbox", {
        plotSpatialFeature(sfe, "nCounts",
                           annotGeometryName = "myofiber_simplified",
                           annot_aes = list(fill = "area"), bbox = bbox_2s1)
    })
    expect_ggplot("Two samples, only colGeometry, different bbox", {
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s)
    })
    expect_ggplot("Two samples, with annotGeometry, different bbox", {
        plotSpatialFeature(sfe, "nCounts",
                           annotGeometryName = "myofiber_simplified",
                           annot_aes = list(fill = "area"), bbox = bbox_2s)
    })
})

colGraphs(sfe, sample_id = "all", name = "visium") <- findVisiumGraph(sfe, sample_id = "all")
sfe <- colDataUnivariate(sfe, "localmoran", features = "nCounts", sample_id = "all")
annotGraphs(sfe, sample = "all", name = "knn") <-
    findSpatialNeighbors(sfe, MARGIN = 3, type = "myofiber_simplified",
                         method = "knearneigh", k = 5, zero.policy = TRUE,
                         sample_id = "all")
sfe <- annotGeometryUnivariate(sfe, "localmoran", "area",
                               annotGeometryName = "myofiber_simplified",
                               zero.policy = TRUE, sample_id = "all")


sfe_muscle <- colDataUnivariate(sfe_muscle, "localmoran", features = "nCounts")
annotGraph(sfe_muscle, "poly") <- findSpatialNeighbors(sfe_muscle, MARGIN = 3,
                                                       type = "myofiber_simplified",
                                                       method = "poly2nb",
                                                       zero.policy = TRUE)
sfe_muscle <- annotGeometryUnivariate(sfe_muscle, "localmoran", "area",
                                      annotGeometryName = "myofiber_simplified",
                                      zero.policy = TRUE)

test_that("Using bbox with plotLocalResults", {
    expect_ggplot("One sample, colGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "nCounts", bbox = bbox2,
                        colGeometryName = "spotPoly", divergent = TRUE,
                        diverge_center = 0)
    })
    expect_ggplot("One sample col and annotGeometry plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "nCounts", bbox = bbox2,
                        colGeometryName = "spotPoly",
                        annotGeometryName = "myofiber_simplified",
                        annot_fixed = list(linewidth = 0.3),
                        divergent = TRUE, diverge_center = 0)
    })
    expect_ggplot("One sample, annotGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "area", bbox = bbox2,
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
    expect_ggplot("Two samples, colGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "nCounts", bbox = bbox_2s1,
                        colGeometryName = "spotPoly",
                        divergent = TRUE, diverge_center = 0)
    })
    expect_ggplot("Two samples col and annotGeometry plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "nCounts", bbox = bbox_2s1,
                        colGeometryName = "spotPoly",
                        annotGeometryName = "myofiber_simplified",
                        annot_fixed = list(linewidth = 0.3),
                        divergent = TRUE, diverge_center = 0)
    })
    expect_ggplot("Two samples, annotGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "area", bbox = bbox_2s1,
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
})
sfe <- logNormCounts(sfe)
sfe <- runPCA(sfe, ncomponents = 20, BSPARAM = BiocSingular::ExactParam())
test_that("Using bbox with spatialReducedDim", {
    expect_ggplot("Use bbox with spatialReducedDim", {
        spatialReducedDim(sfe, dimred = "PCA", ncomponents = 1,
                          colGeometryName = "spotPoly", bbox = bbox_2s1,
                          divergent = TRUE, diverge_center = 0)
    })
    expect_ggplot("Use bbox with spatialReducedDim, 2 PCs", {
        spatialReducedDim(sfe, dimred = "PCA", ncomponents = 2,
                          colGeometryName = "spotPoly", bbox = bbox_2s1,
                          divergent = TRUE, diverge_center = 0, ncol = 1)
    })
})

test_that("Plot graphs with bbox", {
    expect_ggplot("colGraph with bbox", {
        plotColGraph(sfe, colGraphName = "visium",
                     colGeometryName = "spotPoly", bbox = bbox_2s)
    })
    expect_ggplot("annotGraph with bbox", {
        plotAnnotGraph(sfe, annotGraphName = "knn",
                     annotGeometryName = "myofiber_simplified", bbox = bbox_2s)
    })
})

test_that("Using bbox with plotCellBin2D", {
    expect_doppelganger("plotCellBin2D with bbox", {
        plotCellBin2D(sfe_cosmx, bins = 50, bbox = bbox_large)
    })
})

test_that("Incorrect formats of bbox", {
    expect_error(plotSpatialFeature(sfe, "nCounts", bbox = as.list(bbox_2s1)),
                 "numeric vector or a matrix")
    expect_error(plotSpatialFeature(sfe, "nCounts", bbox = "foobar"),
                 "must be a numeric vector")
    expect_error(plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s1[1:3]),
                 "must have length 4")
    expect_error(plotSpatialFeature(sfe, "nCounts",
                                    bbox = setNames(bbox_2s1, c("foo", "bar", "meow", "purr"))),
                 "must be same as the set")
    expect_warning(plotSpatialFeature(sfe, "nCounts",
                                      bbox = unname(bbox_2s1[c("xmin", "ymin", "xmax", "ymax")])),
                   "No names available for bbox. Assuming")
    expect_ggplot("OK if bbox matrix is transposed", {
        plotSpatialFeature(sfe, "nCounts", bbox = t(bbox_2s))
    })
    expect_error(plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s[1:3,]),
                 "must have 4 rows")
    expect_error({
        bbox_test <- bbox_2s
        rownames(bbox_test) <- c("foo", "bar", "meow", "purr")
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_test)
    }, "Row names of bbox must be same as the set")
    expect_error({
        bbox_test <- bbox_2s
        colnames(bbox_test) <- c("foo", "sample02")
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_test)
    }, "Column names of bbox must match the sample IDs")
    expect_error(plotSpatialFeature(sfe, "nCounts", bbox = bbox),
                 "The bounding box does not overlap with the geometry")
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
    expect_doppelganger("Moran plot don't plot influential", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, bins = 30,
                  plot_influential = FALSE)
    })
    expect_doppelganger("Moran plot hex bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, hex = TRUE, bins = 30)
    })
})

test_that("Plot geometries", {
    expect_ggplot("plot colGeometry 2 samples", {
        plotGeometry(sfe, "spotPoly")
    })
    expect_ggplot("plot colGeometry 1 sample", {
        plotGeometry(sfe_muscle, "spotPoly")
    })
    expect_ggplot("plot annotGeometry 2 samples", {
        plotGeometry(sfe, "myofiber_simplified", MARGIN = 3)
    })
    expect_ggplot("plot annotGeometry 1 sample", {
        plotGeometry(sfe_muscle, "myofiber_simplified", MARGIN = 3)
    })
    expect_ggplot("Plot colGeometry, with bbox", {
        plotGeometry(sfe, "spotPoly", bbox = bbox_2s)
    })
    expect_ggplot("Plot annotGeometry, with bbox", {
        plotGeometry(sfe, "myofiber_simplified", MARGIN = 3, bbox = bbox_2s)
    })
})

test_that("Message about using linewidth instead of size for polygon outlines", {
    expect_message(plotSpatialFeature(sfe, "nCounts", fill = NA, size = 0.5,
                                      aes_use = "color"),
                   "Please use linewidth instead of size for thickness of polygon outlines.")
    # Still get the right plot
    expect_ggplot("Plot polygon, with size rather than linewidth",
                        plotSpatialFeature(sfe, "nCounts", fill = NA, size = 0.5,
                                           aes_use = "color"))
    expect_doppelganger("Moran plot hex bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, hex = TRUE, bins = 30)
    })
})

sfe_muscle2 <- logNormCounts(sfe_muscle2)
gs <- modelGeneVar(sfe_muscle2)
hvgs <- getTopHVGs(gs, fdr.threshold = 0.05)
sfe_muscle2 <- runMultivariate(sfe_muscle2, "multispati", subset_row = hvgs)
sfe_muscle2 <- reducedDimUnivariate(sfe_muscle2, "sp.correlogram", dimred = "multispati",
                                    components = 1:10, order = 3)
set.seed(29)
sfe_muscle2 <- reducedDimUnivariate(sfe_muscle2, "moran.mc", dimred = "multispati",
                                    components = 1:5, nsim = 99)
sfe_muscle2 <- reducedDimUnivariate(sfe_muscle2, "moran.plot", dimred = "multispati",
                                   components = 1, colGraphName = "visium")
test_that("Univariate downstream plots for dimred", {
    expect_ggplot("Moran plot for PCA", {
        suppressWarnings(moranPlot(sfe_muscle2, "PC1", graphName = "visium",
                                   reducedDimName = "PCA"))
    })
    expect_doppelganger("Correlograms for multispati PCs", {
        plotCorrelogram(sfe_muscle2, 1:10, reducedDimName = "multispati")
    })
    expect_doppelganger("Correlograms for multispati PCs, one component", {
        plotCorrelogram(sfe_muscle2, "PC1", reducedDimName = "multispati")
    })
})

test_that("Moran MC plot for dimred", {
    expect_ggplot("Use dim names", {
        plotMoranMC(sfe_muscle2, "PC1", reducedDimName = "multispati")
    })
    expect_ggplot("Use numeric indices", {
        plotMoranMC(sfe_muscle2, 1:5, reducedDimName = "multispati")
    })
})

# Plot image behind geometries
# Need uncropped image
if (!dir.exists("outs")) dir.create("outs")
mat_fn <- file.path("outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
                  destfile = file.path("outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
                  destfile = file.path("outs", "spatial.tar.gz"))
    untar(file.path("outs", "spatial.tar.gz"), exdir = "outs")
}

sfe_ob <- read10xVisiumSFE(".", images = "lowres")

# get 2 samples
inds <- spatialCoords(sfe_ob)[,1] < 5000
sfe_ob1 <- sfe_ob[,inds]
sfe_ob2 <- sfe_ob[,!inds]
sfe_ob2 <- changeSampleIDs(sfe_ob2, c(sample01 = "sample02"))
sfe_ob3 <- SpatialFeatureExperiment::cbind(sfe_ob1, sfe_ob2)
sfe_ob3 <- removeEmptySpace(sfe_ob3)

bbox_use <- c(xmin = 4000, ymin = 6000, xmax = 4750, ymax = 6750)

dir_use <- system.file(file.path("extdata", "vizgen"), package = "SpatialFeatureExperiment")
sfe_mer <- readVizgen(dir_use, z = 0L, image = "PolyT", use_cellpose = FALSE)

test_that("plotSpatialFeature with RGB image in the background", {
    expect_ggplot("One sample, one feature", {
        plotSpatialFeature(sfe_ob, "array_row", image_id = "lowres")
    })
    expect_ggplot("One sample, two features", {
        plotSpatialFeature(sfe_ob, c("in_tissue", "array_col"),
                           image_id = "lowres", maxcell = 5e+4)
    })
    expect_ggplot("With bbox", {
        plotSpatialFeature(sfe_ob, "array_row", image_id = "lowres",
                           bbox = bbox_use)
    })
    expect_ggplot("Two samples, one feature", {
        plotSpatialFeature(sfe_ob3, "array_row", image_id = "lowres")
    })
    expect_ggplot("Two samples, two features", {
        plotSpatialFeature(sfe_ob3, c("in_tissue", "array_col"),
                           image_id = "lowres", maxcell = 5e+4)
    })
    expect_ggplot("One sample, one feature, grayscale", {
        plotSpatialFeature(sfe_mer, "volume", image_id = "PolyT",
                           colGeometryName = "cellSeg", alpha = 0.5,
                           dark = TRUE)
    })
})

colGraph(sfe_ob, "visium") <- findVisiumGraph(sfe_ob)
sfe_ob <- logNormCounts(sfe_ob)
gs <- modelGeneVar(sfe_ob)
hvgs <- getTopHVGs(gs, fdr.threshold = 0.01)
sfe_ob <- runUnivariate(sfe_ob, "localmoran", hvgs[1])

ag <- spotPoly(sfe_ob)
ag$gene <- logcounts(sfe_ob)[hvgs[2],]
annotGeometry(sfe_ob, "foo") <- ag
annotGraph(sfe_ob, "bar") <- colGraph(sfe_ob, "visium")
sfe_ob <- annotGeometryUnivariate(sfe_ob, "LOSH", "gene", annotGeometryName = "foo",
                                  annotGraphName = "bar")

test_that("plotLocalResult with image", {
    expect_ggplot("colGeometry", {
        plotLocalResult(sfe_ob, "localmoran", hvgs[1], swap_rownames = "symbol",
                        colGeometryName = "spotPoly", divergent = TRUE,
                        diverge_center = 0, image_id = "lowres")
    })
    expect_ggplot("annotGeometry", {
        plotLocalResult(sfe_ob, "LOSH", "gene",
                        annotGeometryName = "foo", image_id = "lowres")
    })
})

sfe_ob <- runPCA(sfe_ob, ncomponents = 5, subset_row = hvgs)
test_that("spatialReducedDim with image", {
    expect_ggplot("", {
        spatialReducedDim(sfe_ob, "PCA", 2, image_id = "lowres", maxcell = 5e4,
                          divergent = TRUE, diverge_center = 0)
    })
})

test_that("plotGeometry with image", {
    expect_ggplot("One sample", {
        plotGeometry(sfe_ob, "spotPoly", MARGIN = 2, image_id = "lowres")
    })
    expect_ggplot("Two samples", {
        plotGeometry(sfe_ob3, "spotPoly", MARGIN = 2, image_id = "lowres")
    })
})
