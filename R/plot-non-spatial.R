#' Plot the elbow plot or scree plot for PCA
#'
#' Apparently, there is no apparent way to plot the PC elbow plot other than
#' extracting the variance explained attribute of the dimred slot, because even
#' the OSCA book makes the elbow plot this way, which I find kind of cumbersome
#' compared to Seurat. So I'm writing this function to make the elbow plot with
#' SCE less cumbersome.
#'
#' @param sce A \code{SingleCellExperiment} object, or anything that inherits
#' from \code{SingleCellExperiment}.
#' @param ndims Number of PCs to plot.
#' @param reduction Name of the dimension reduction to use. It must have an
#' attribute called "percentVar". Defaults to "PCA".
#' @return A ggplot object. The y axis is percentage of variance explained.
#' @importFrom scales breaks_extended
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- runPCA(sfe, ncomponents = 10, exprs_values = "counts")
#' ElbowPlot(sfe, ndims = 10)
ElbowPlot <- function(sce, ndims = 20, reduction = "PCA") {
    # For scater::runPCA results
    percent_var <- attr(reducedDim(sce, reduction), "percentVar")
    if (length(percent_var) < ndims) ndims <- length(percent_var)
    inds <- seq_len(ndims)
    df <- data.frame(
        PC = inds,
        pct_var = percent_var[inds]
    )
    PC <- pct_var <- NULL
    ggplot(df, aes(PC, pct_var)) +
        geom_point() +
        labs(x = "PC", y = "Variance explained (%)") +
        scale_x_continuous(breaks = breaks_extended(n = min(ndims, 10), Q = 1:5)) +
        theme(panel.grid.minor.x = element_blank())
}

.get_top_loading_genes <- function(df, nfeatures, balanced) {
    # df has columns gene_show and value
    if (balanced) {
        n2 <- floor(nfeatures / 2)
        ord_plus <- order(df$value, decreasing = TRUE)
        ord_minus <- order(df$value, decreasing = FALSE)
        out <- rbind(df[ord_plus[seq_len(n2)], ], df[ord_minus[seq_len(n2)], ])
    } else {
        ord <- order(abs(df$value), decreasing = TRUE)
        out <- df[ord[seq_len(nfeatures)], ]
    }
    return(out)
}

#' Plot top PC loadings of genes
#'
#' Just like Seurat's VizDimLoadings function. I haven't found an equivalent for
#' SCE but find it useful. But I'm not trying to reproduce that Seurat function
#' exactly. For instance, I don't like it when Seurat imposes a ggplot theme,
#' and I don't like the cowplot theme. Maybe I should rewrite it in base R but
#' for now I'm using Tidyverse.
#'
#' @inheritParams ElbowPlot
#' @param dims Numeric vector specifying which PCs to plot.
#' @param nfeatures Number of genes to plot.
#' @param show_symbol Logical; if the row names of the matrix are Ensembl
#'   accessions, indicate whether to show more human readable gene symbols in
#'   the plot instead. Ignored if the column specified in \code{symbol_col} is
#'   absent from rowData.
#' @param symbol_col If the row names of the gene expression matrix are Ensembl
#'   accessions to avoid ambiguity in analysis. If not found in \code{rowData},
#'   then rownames of the gene count matrix will be used.
#' @param balanced Return an equal number of genes with + and - scores. If
#'   FALSE, returns the top genes ranked by the scores absolute values.
#' @param ncol Number of columns in the facetted plot.
#' @return A ggplot object. Loadings for different PCs are plotted in different
#'   facets so one ggplot object is returned.
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 facet_wrap scale_y_discrete
#' @importFrom stats reorder
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- runPCA(sfe, ncomponents = 10, exprs_values = "counts")
#' plotDimLoadings(sfe, dims = 1:2)
plotDimLoadings <- function(sce, dims = 1:4, nfeatures = 10,
                            show_symbol = TRUE, symbol_col = "symbol",
                            reduction = "PCA",
                            balanced = TRUE, ncol = 2) {
    # For scater::runPCA results
    loadings <- attr(reducedDim(sce, reduction), "rotation")
    is_ensembl <- all(grepl("^ENS", rownames(sce)))
    if (!is_ensembl) show_symbol <- FALSE # I mean, irrelevant
    loading_cols <- paste0("PC", dims)
    df <- cbind(as.data.frame(rowData(sce)[rownames(loadings),, drop = FALSE]),
                loadings[, loading_cols])
    if (!symbol_col %in% names(df) || !show_symbol) {
        df$gene_show <- rownames(loadings)
    } else {
        df$gene_show <- df[[symbol_col]]
    }
    df_plt <- lapply(loading_cols, function(p) {
        df_use <- df[, c("gene_show", p)]
        names(df_use)[2] <- "value"
        out <- .get_top_loading_genes(df_use, nfeatures, balanced)
        out$PC <- p
        out
    })
    df_plt <- Reduce(rbind, df_plt)
    df_plt$PC <- factor(df_plt$PC, levels = loading_cols)
    # Basically reimplementing tidytext::reorder_within and scale_y_reordered
    df_plt$gene <- paste(df_plt$gene_show, df_plt$PC, sep = "___")
    df_plt$gene <- reorder(df_plt$gene, df_plt$value)
    reg <- "___.+$"
    value <- gene <- NULL
    ggplot(df_plt, aes(value, gene)) +
        geom_segment(aes(yend = gene), xend = 0, show.legend = FALSE) +
        geom_point(color = "blue") +
        geom_vline(xintercept = 0, linetype = 2) +
        facet_wrap(~PC, scales = "free_y", ncol = ncol) +
        scale_y_discrete(labels = function(x) gsub(reg, "", x)) +
        labs(x = "Loading", y = "Gene")
}

.plot_dimdata_bin2d_fun <- function(fun) {
    function(sfe, x, y, subset = NULL, bins = 100, binwidth = NULL,
             hex = FALSE, name_true = NULL, name_false = NULL) {
        bin_fun <- if (hex) geom_hex else geom_bin2d
        df <- as.data.frame(fun(sfe))
        p <- ggplot()
        if (is.null(subset)) {
            p <- p +
                bin_fun(aes(.data[[x]], .data[[y]]), bins = bins, binwidth = binwidth,
                        data = df) +
                scale_fill_distiller(palette = "Blues", direction = 1)
        } else {
            name_subset <- subset
            subset <- df[[subset]]
            if (anyNA(as.logical(subset))) {
                stop("Column ", subset, " must be coerceable to logical.")
            }
            name_true <- name_true %||% name_subset
            name_false <- name_false %||% paste0("!", name_subset)
            p <- p +
                bin_fun(aes(.data[[x]], .data[[y]]), bins = bins, binwidth = binwidth,
                        data = df[!subset,]) +
                scale_fill_distiller(palette = "Blues", direction = 1,
                                     name = name_false) +
                new_scale_fill() +
                bin_fun(aes(.data[[x]], .data[[y]]), bins = bins, binwidth = binwidth,
                        data = df[subset,]) +
                scale_fill_distiller(palette = "Reds", direction = 1,
                                     name = name_true)
        }
        p
    }
}

#' Plot colData and rowData with 2D histograms
#'
#' To avoid overplotting in large datasets. The 2D histogram is more informative
#' of point density on the plot than the scatter plot where there are so many
#' points plotted that they effectively form a solid block.
#'
#' @inheritParams ggplot2::geom_bin2d
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param x Name of the column in \code{colData} or \code{rowData} to plot on
#'   the x axis of the plot.
#' @param y Name of the column in \code{colData} or \code{rowData} to plot on
#'   the y axis of the plot.
#' @param bins Numeric vector giving number of bins in both vertical and
#'   horizontal directions. Set to 100 by default.
#' @param subset Name of a logical column in \code{colData} or \code{rowData},
#'   indicating cells or genes to plot with a different palette. Since the 2D
#'   histogram is effectively an opaque heatmap, don't use this argument unless
#'   the two groups are largely non-overlapping in the variables being plotted.
#' @param hex Logical, whether to use hexagon rather than rectangular bins.
#'   Requires the \code{hexbin} package.
#' @param name_true Character, name to show on the legend for cells or genes
#'   indicated \code{TRUE} in the \code{subset} argument.
#' @param name_false Character, name to show on the legend for cells or genes
#'   indicated \code{FALSE} in the \code{subset} argument.
#' @importFrom ggplot2 geom_bin2d geom_hex
#' @importFrom stats reshape
#' @export
#' @return A ggplot object
#' @name plotColDataBin2D
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData()
#' sfe <- sfe[, sfe$in_tissue]
#' plotColDataBin2D(sfe, "nCounts", "nGenes")
plotColDataBin2D <- .plot_dimdata_bin2d_fun(colData)

#' @rdname plotColDataBin2D
#' @export
plotRowDataBin2D <- .plot_dimdata_bin2d_fun(rowData)

.plot_dimdata_hist <- function(fun) {
    function(sfe, feature, fill_by = NULL, subset = NULL, bins = 100, binwidth = NULL,
             scales = "free", ncol = 1, position = "identity") {
        df <- as.data.frame(fun(sfe))[, c(feature, fill_by, subset), drop = FALSE]
        if (!is.null(subset)) df <- df[df[[subset]],]
        p <- ggplot()
        if (length(feature) > 1L) {
            df <- reshape(df, varying = feature, direction = "long",
                          v.names = "values", timevar = "variable",
                          times = feature)
            p <- p +
                geom_histogram(data = df,
                               mapping = aes(!!!syms(c(x = "values", fill = fill_by))),
                               bins = bins,
                               binwidth = binwidth, position = position) +
                facet_wrap(~ variable, scales = scales, ncol = ncol)
        } else {
            p <- p +
                geom_histogram(data = df,
                               mapping = aes(!!!syms(c(x = feature, fill = fill_by))),
                               bins = bins,
                               binwidth = binwidth, position = position)
        }
        p
    }
}

#' Plot histograms for colData and rowData columns
#'
#' @inheritParams plotColDataBin2D
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams ggplot2::geom_histogram
#' @param feature Names of columns in \code{colData} or \code{rowData} to plot.
#'   When multiple features are specified, they will be plotted in separate
#'   facets.
#' @param ncol Number of columns in the facetting.
#' @param fill_by Name of a categorical column in \code{colData} or
#'   \code{rowData} to fill the histogram.
#' @param subset Name of a logical column to only plot a subset of the data.
#' @return A ggplot object
#' @seealso plotColDataFreqpoly
#' @importFrom rlang %||% .data
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData()
#' plotColDataHistogram(sfe, c("nCounts", "nGenes"), fill_by = "in_tissue",
#'                      bins = 50, position = "stack")
#' plotColDataHistogram(sfe, "nCounts", subset = "in_tissue")
#' sfe2 <- sfe[, sfe$in_tissue]
#' plotColDataHistogram(sfe2, c("nCounts", "nGenes"), bins = 50)
plotColDataHistogram <- .plot_dimdata_hist(colData)

#' @rdname plotColDataHistogram
#' @export
plotRowDataHistogram <- .plot_dimdata_hist(rowData)

.plot_dimdata_freqpoly <- function(fun) {
    function(sfe, feature, color_by = NULL, subset = NULL, bins = 100,
             binwidth = NULL, linewidth = 1.2,
             scales = "free", ncol = 1, position = "identity") {
        df <- as.data.frame(fun(sfe))[, c(feature, color_by, subset), drop = FALSE]
        if (!is.null(subset)) df <- df[df[[subset]],]
        p <- ggplot()
        if (length(feature) > 1L) {
            df <- reshape(df, varying = feature, direction = "long",
                          v.names = "values", timevar = "variable",
                          times = feature)
            p <- p +
                geom_freqpoly(data = df,
                              mapping = aes(!!!syms(c(x = "values", color = color_by))),
                              bins = bins, linewidth = linewidth,
                              binwidth = binwidth, position = position) +
                facet_wrap(~ variable, scales = scales, ncol = ncol)
        } else {
            p <- p +
                geom_freqpoly(data = df,
                              mapping = aes(!!!syms(c(x = feature, color = color_by))),
                              bins = bins, linewidth = linewidth,
                              binwidth = binwidth, position = position)
        }
        p <- p + scale_color_manual(values = ditto_colors)
        p
    }
}

#' Plot frequency polygons for colData and rowData columns
#'
#' This function is recommended instead of \code{\link{plotColDataHistogram}}
#' when coloring by multiple categories and log transforming the y axis, which
#' causes problems in stacked histograms.
#'
#' @inheritParams ggplot2::geom_freqpoly
#' @inheritParams plotColDataHistogram
#' @param linewidth Line width of the polygons, defaults to a thicker 1.2.
#' @param color_by Name of a categorical column in \code{colData} or
#'   \code{rowData} to color the polygons.
#' @seealso plotColDataHistogram
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData()
#' plotColDataFreqpoly(sfe, c("nCounts", "nGenes"), color_by = "in_tissue",
#'                     bins = 50)
#' plotColDataFreqpoly(sfe, "nCounts", subset = "in_tissue")
#' sfe2 <- sfe[, sfe$in_tissue]
#' plotColDataFreqpoly(sfe2, c("nCounts", "nGenes"), bins = 50)
plotColDataFreqpoly <- .plot_dimdata_freqpoly(colData)

#' @rdname plotColDataFreqpoly
#' @export
plotRowDataFreqpoly <- .plot_dimdata_freqpoly(rowData)
