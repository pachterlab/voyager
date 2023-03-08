# Unit test plotting functions with vdiffr
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(vdiffr)
library(scater)
library(Matrix)
library(ggplot2)

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
    expect_doppelganger("Plot localmoran Ii for gene", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_doppelganger("Plot Ii for gene on top of an annotation", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            colGeometryName = "spotPoly",
            annotGeometryName = "tissueBoundary", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_doppelganger("Plot another column", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "Z.Ii",
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_doppelganger("Plot a categorical attribute", {
        plotLocalResult(sfe_muscle, "localmoran", "Myh1",
            attribute = "mean",
            colGeometryName = "spotPoly", swap_rownames = "symbol"
        )
    })
    expect_doppelganger("Plot 2 features", {
        plotLocalResult(sfe_muscle, "localmoran", c("Myh1", "Myh2"),
            colGeometryName = "spotPoly", divergent = TRUE,
            diverge_center = 0, swap_rownames = "symbol"
        )
    })
    expect_doppelganger("Plot Ii for annotGeometry alone", {
        plotLocalResult(sfe_muscle, "localmoran", "area", "Ii",
            annotGeometryName = "myofiber_simplified",
            linewidth = 0.3, color = "cyan", divergent = TRUE,
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
    skip_on_ci()
    expect_warning(
        moranPlot(sfe, "B", "visium", "sample01"),
        "Too few points"
    )
    expect_doppelganger(
        "moranPlot, not filled",
        moranPlot(sfe_muscle, feature_use, "visium", swap_rownames = "symbol")
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
            color_by = "GraphBased", contour_color = "blue",
            swap_rownames = "symbol"
        )
    )
})

test_that("moranPlot, filled, no color_by", {
    skip_on_ci()
    expect_doppelganger(
        "moranPlot, filled",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE, swap_rownames = "symbol"
        )
    )
})

test_that("moranPlot, filled, with color_by", {
    skip_on_ci()
    expect_doppelganger(
        "moranPlot, filled, with color_by",
        moranPlot(sfe_muscle, feature_use, "visium",
            filled = TRUE, color_by = "GraphBased", swap_rownames = "symbol"
        )
    )
})

test_that("plot graphs", {
    expect_doppelganger(
        "plotColGraph",
        plotColGraph(sfe_muscle, colGraphName = "visium",
                     colGeometryName = "spotPoly")
    )
    expect_doppelganger("plotAnnotGraph",
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
                        plotSpatialFeature(sfe_muscle, "HLA-foo", "spotPoly"))
    sfe_muscle <- runUnivariate(sfe_muscle, "localmoran", "HLA-foo",
                                colGraphName = "visium")
    expect_doppelganger("plotLocalResult with illegal gene name",
                        plotLocalResult(sfe_muscle, "localmoran", "HLA-foo",
                                        colGeometryName = "spotPoly"))
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
    expect_message(plotSpatialFeature(sfe_cosmx, "nCounts",
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
    expect_warning(plotColDataBin2D(sfe = sfe_cosmx, x = "nCounts",
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
})

test_that("plotCellBin2D", {
    expect_doppelganger("Cell density, rectangular", {
        plotCellBin2D(sfe_cosmx, bins = 50)
    })
<<<<<<< HEAD
    expect_doppelganger("Cell density, hex", {
        plotCellBin2D(sfe_cosmx, hex = TRUE, bins = 50)
    })
    expect_doppelganger("Multiple samples", {
        plotCellBin2D(sfe_cosmx2, bins = 50)
    })
})

test_that("Binning values", {
    expect_doppelganger("Bin and summarize feature", {
        plotSpatialFeature(sfe_cosmx, features = "nCounts",
                           colGeometryName = "centroids", bins = 50)
    })
    expect_doppelganger("Bin and summarize dimension reduction values", {
        spatialReducedDim(sfe_cosmx, dimred = "PCA", ncomponents = 1,
                          colGeometryName = "centroids", bins = 50)
    })
    expect_doppelganger("Bin and summarize local results", {
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
    expect_doppelganger("Only plotting colGeometry", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "cellSeg",
                           bbox = bbox)
    })
    expect_doppelganger("With scattermore", {
        plotSpatialFeature(sfe_cosmx, "nCounts", colGeometryName = "centroids",
                           bbox = bbox_large, scattermore = TRUE, pointsize = 1)
    })
    expect_doppelganger("Both colGeometry and annotGeometry", {
        plotSpatialFeature(sfe_muscle, "nCounts",
                           annotGeometryName = "myofiber_simplified",
                           annot_aes = list(fill = "area"), bbox = bbox2)
    })
    # Two samples
    expect_doppelganger("Two samples, only plotting colGeometry, same bbox", {
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s1)
    })
    expect_doppelganger("Two samples, with annotGeometry, same bbox", {
        plotSpatialFeature(sfe, "nCounts",
                           annotGeometryName = "myofiber_simplified",
                           annot_aes = list(fill = "area"), bbox = bbox_2s1)
    })
    expect_doppelganger("Two samples, only colGeometry, different bbox", {
        plotSpatialFeature(sfe, "nCounts", bbox = bbox_2s)
    })
    expect_doppelganger("Two samples, with annotGeometry, different bbox", {
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
    expect_doppelganger("One sample, colGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "nCounts", bbox = bbox2,
                        colGeometryName = "spotPoly", divergent = TRUE,
                        diverge_center = 0)
    })
    expect_doppelganger("One sample col and annotGeometry plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "nCounts", bbox = bbox2,
                        colGeometryName = "spotPoly",
                        annotGeometryName = "myofiber_simplified",
                        annot_fixed = list(linewidth = 0.3),
                        divergent = TRUE, diverge_center = 0)
    })
    expect_doppelganger("One sample, annotGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe_muscle, "localmoran", "area", bbox = bbox2,
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
    expect_doppelganger("Two samples, colGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "nCounts", bbox = bbox_2s1,
                        colGeometryName = "spotPoly",
                        divergent = TRUE, diverge_center = 0)
    })
    expect_doppelganger("Two samples col and annotGeometry plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "nCounts", bbox = bbox_2s1,
                        colGeometryName = "spotPoly",
                        annotGeometryName = "myofiber_simplified",
                        annot_fixed = list(linewidth = 0.3),
                        divergent = TRUE, diverge_center = 0)
    })
    expect_doppelganger("Two samples, annotGeometry, plotLocalResults bbox", {
        plotLocalResult(sfe, "localmoran", "area", bbox = bbox_2s1,
                        annotGeometryName = "myofiber_simplified",
                        divergent = TRUE, diverge_center = 0)
    })
})
sfe <- logNormCounts(sfe)
sfe <- runPCA(sfe, ncomponents = 20, BSPARAM = BiocSingular::ExactParam())
test_that("Using bbox with spatialReducedDim", {
    expect_doppelganger("Use bbox with spatialReducedDim", {
        spatialReducedDim(sfe, dimred = "PCA", ncomponents = 1,
                          colGeometryName = "spotPoly", bbox = bbox_2s1,
                          divergent = TRUE, diverge_center = 0)
    })
    expect_doppelganger("Use bbox with spatialReducedDim, 2 PCs", {
        spatialReducedDim(sfe, dimred = "PCA", ncomponents = 2,
                          colGeometryName = "spotPoly", bbox = bbox_2s1,
                          divergent = TRUE, diverge_center = 0, ncol = 1)
    })
})

test_that("Plot graphs with bbox", {
    expect_doppelganger("colGraph with bbox", {
        plotColGraph(sfe, colGraphName = "visium",
                     colGeometryName = "spotPoly", bbox = bbox_2s)
    })
    expect_doppelganger("annotGraph with bbox", {
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
    expect_doppelganger("OK if bbox matrix is transposed", {
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
<<<<<<< HEAD
    expect_doppelganger("Moran plot hex bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, hex = TRUE, bins = 30)
    })
})

test_that("Plot geometries", {
    expect_doppelganger("plot colGeometry 2 samples", {
        plotGeometry(sfe, "spotPoly")
    })
    expect_doppelganger("plot colGeometry 1 sample", {
        plotGeometry(sfe_muscle, "spotPoly")
    })
    expect_doppelganger("plot annotGeometry 2 samples", {
        plotGeometry(sfe, "myofiber_simplified", MARGIN = 3)
    })
    expect_doppelganger("plot annotGeometry 1 sample", {
        plotGeometry(sfe_muscle, "myofiber_simplified", MARGIN = 3)
    })
    expect_doppelganger("Plot colGeometry, with bbox", {
        plotGeometry(sfe, "spotPoly", bbox = bbox_2s)
    })
    expect_doppelganger("Plot annotGeometry, with bbox", {
        plotGeometry(sfe, "myofiber_simplified", MARGIN = 3, bbox = bbox_2s)
    })
})

test_that("Message about using linewidth instead of size for polygon outlines", {
    expect_message(plotSpatialFeature(sfe, "nCounts", fill = NA, size = 0.5,
                                      aes_use = "color"),
                   "Please use linewidth instead of size for thickness of polygon outlines.")
    # Still get the right plot
    expect_doppelganger("Plot polygon, with size rather than linewidth",
                        plotSpatialFeature(sfe, "nCounts", fill = NA, size = 0.5,
                                           aes_use = "color"))
    expect_doppelganger("Moran plot hex bin", {
        moranPlot(sfe_muscle2, "nCounts", binned = TRUE, hex = TRUE, bins = 30)
    })
})
