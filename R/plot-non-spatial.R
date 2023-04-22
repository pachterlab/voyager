#' Plot the elbow plot or scree plot for PCA
#'
#' Apparently, there is no apparent way to plot the PC elbow plot other than
#' extracting the variance explained attribute of the dimred slot, because even
#' the OSCA book makes the elbow plot this way, which I find kind of cumbersome
#' compared to Seurat. So I'm writing this function to make the elbow plot with
#' SCE less cumbersome.
#'
#' @inheritParams calculateUnivariate
#' @inheritParams multispati_rsp
#' @param sce A \code{SingleCellExperiment} object, or anything that inherits
#'   from \code{SingleCellExperiment}.
#' @param ndims Number of components with positive eigenvalues, such as PCs in
#'   non-spatial PCA.
#' @param reduction Name of the dimension reduction to use. It must have an
#'   attribute called either "percentVar" or "eig" for eigenvalues. Defaults to
#'   "PCA".
#' @param facet Logical, whether to facet by samples when multiple samples are
#'   present. This is relevant when spatial PCA is run separately for each
#'   sample, which gives different results from running jointly for all samples.
#' @param ncol Number of columns of facets if facetting.
#' @return A ggplot object. The y axis is eigenvalues or percentage of variance
#'   explained if relevant.
#' @importFrom scales breaks_extended
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- runPCA(sfe, ncomponents = 10, exprs_values = "counts")
#' ElbowPlot(sfe, ndims = 10)
ElbowPlot <- function(sce, ndims = 20, nfnega = 0, reduction = "PCA",
                      sample_id = "all", facet = FALSE, ncol = NULL) {
    # For scater::runPCA results
    # to do: 1. deal with other dimension reductions with eigenvalues
    # and negative eigenvalues
    # 2. deal with multiple samples with separate spatial dimred results
    dimred <- reducedDim(sce, reduction)
    use_pct_var <- "percentVar" %in% names(attributes(dimred))
    if (use_pct_var)
        y <- attr(dimred, "percentVar") / 100
    else y <- attr(dimred, "eig")
    if (is.vector(y)) y <- matrix(y, ncol = 1, dimnames = list(NULL, "value"))
    if (sample_id == "all") sample_id <- colnames(y)
    y <- y[,sample_id, drop = FALSE]
    nf <- nrow(y)
    eigs_sign <- rowSums(y)
    ndims <- min(ndims, sum(eigs_sign > 0))
    nfnega <- min(nfnega, sum(eigs_sign < 0))
    inds_posi <- seq_len(ndims)
    inds_nega <- tail(seq_len(nrow(y)), nfnega)
    labels <- as.character(c(inds_posi, inds_nega))
    y_posi <- y[inds_posi,,drop = FALSE]
    y_nega <- y[inds_nega,,drop = FALSE]
    y <- rbind(y_posi, y_nega)
    inds <- seq_len(ndims + nfnega)
    df <- data.frame(PC = inds)
    df <- cbind(df, y)
    if (ncol(y) > 1L) {
        df <- reshape(df, varying = setdiff(colnames(y), "sign"), direction = "long",
                      v.names = "value", timevar = "sample",
                      times = colnames(y))
    }
    df$sign <- ifelse(df$PC > ndims, "n", "p")
    PC <- pct_var <- value <- sign <- NULL
    aes_use <- aes(PC, value)
    use_break <- length(unique(df$sign)) > 1L
    do_color <- !(facet || ncol(y) == 1L)
    if (do_color) {
        aes_use <- modifyList(aes_use, aes(color = sample))
    } else if (use_break)
        aes_use <- modifyList(aes_use, aes(group = sign))
    p <- ggplot(df, aes_use)
    if (do_color)
        p <- p + scale_color_manual(values = ditto_colors)

    breaks_inds <- breaks_extended(n = min(ndims+nfnega, 10), Q = 1:5)(inds)
    labels_use <- labels[breaks_inds]
    if (breaks_inds[1] == 0) labels_use <- c("", labels_use)
    if (do_color && use_break) {
        # Plot the positive and negative parts with separate data frames
        df1 <- df[df$sign == "p",]
        df2 <- df[df$sign == "n",]
        p <- p + geom_line(data = df1) + geom_line(data = df2)
    } else p <- p + geom_line()
    p <- p +
        geom_point() +
        scale_x_continuous(breaks = breaks_inds, labels = labels_use) +
        theme(panel.grid.minor.x = element_blank())
    if (nfnega > 0 && ndims > 0) {
        p <- p + geom_hline(yintercept = 0, linetype = 2)
    }
    if (use_pct_var) {
        p <- p + labs(x = "PC", y = "Variance explained") +
            scale_y_continuous(labels = scales::percent)
    } else {
        p <- p + labs(x = "PC", y = "Eigenvalue")
    }
    if (facet) p <- p + facet_wrap(~ sample, ncol = ncol)
    p
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
# For each sample
.get_loadings_df <- function(sce, loadings, loading_cols, nfeatures, balanced,
                             swap_rownames) {
    df <- cbind(as.data.frame(rowData(sce)[rownames(loadings),, drop = FALSE]),
                loadings[, loading_cols])
    if (is.null(swap_rownames) || !swap_rownames %in% names(df)) {
        df$gene_show <- rownames(loadings)
    } else {
        df$gene_show <- df[[swap_rownames]]
    }
    df_plt <- lapply(loading_cols, function(p) {
        df_use <- df[, c("gene_show", p)]
        names(df_use)[2] <- "value"
        out <- .get_top_loading_genes(df_use, nfeatures, balanced)
        out$PC <- p
        out
    })
    df_plt <- Reduce(rbind, df_plt)
    df_plt
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
#' @param dims Numeric vector specifying which PCs to plot. For MULTISPATI, PCs
#' with negative eigenvalues are in the right most columns of the embedding and
#' loading matrices. See the \code{\link{ElbowPlot}}.
#' @param nfeatures Number of genes to plot.
#' @param show_symbol Deprecated. Use argument \code{swap_rownames} instead, to
#'   be consistent with \code{scater} plotting functions.
#' @param symbol_col Deprecated. Use argument \code{swap_rownames} instead, to
#'   be consistent with \code{scater} plotting functions.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to
#'   identify features instead of \code{rownames(object)} when labeling plot
#'   elements. If not found in \code{rowData}, then rownames of the gene count
#'   matrix will be used.
#' @param balanced Return an equal number of genes with + and - scores. If
#'   FALSE, returns the top genes ranked by the scores absolute values.
#' @param ncol Number of columns in the facetted plot.
#' @return A ggplot object. Loadings for different PCs are plotted in different
#'   facets so one ggplot object is returned.
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 facet_wrap scale_y_discrete
#' @importFrom stats reorder
#' @importFrom lifecycle deprecated is_present deprecate_warn
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- runPCA(sfe, ncomponents = 10, exprs_values = "counts")
#' plotDimLoadings(sfe, dims = 1:2)
plotDimLoadings <- function(sce, dims = 1:4, nfeatures = 10,
                            swap_rownames = NULL,
                            show_symbol = deprecated(),
                            symbol_col = deprecated(),
                            reduction = "PCA",
                            balanced = TRUE, ncol = 2,
                            sample_id = "all") {
    # deal with multiple samples with separate spatial dimred results
    if (is_present(show_symbol)) {
        deprecate_warn("1.2.0", "plotDimLoadings(show_symbol = )",
                       "plotDimLoadings(swap_rownames = )")
        if (is_present(symbol_col) && is.null(swap_rownames))
            swap_rownames <- symbol_col
    }

    loadings <- attr(reducedDim(sce, reduction), "rotation")
    loading_cols <- paste0("PC", dims)
    is_multi <- length(dim(loadings)) > 2L
    if (is_multi) {
        df_plts <- lapply(dimnames(loadings)[[3]], function(s) {
            o <- .get_loadings_df(sce, loadings[,,s], loading_cols, nfeatures, balanced,
                                  swap_rownames)
            o$sample <- s
            o
        })
        df_plt <- do.call(rbind, df_plts)
    } else {
        df_plt <- .get_loadings_df(sce, loadings, loading_cols, nfeatures, balanced,
                                   swap_rownames)
    }
    df_plt$PC <- factor(df_plt$PC, levels = loading_cols)
    # Basically reimplementing tidytext::reorder_within and scale_y_reordered
    df_plt$gene <- paste(df_plt$gene_show, df_plt$PC, sep = "___")
    if (is_multi)
        df_plt$gene <- paste(df_plt$gene, df_plt$sample, sep = "___")
    df_plt$gene <- reorder(df_plt$gene, df_plt$value)
    reg <- "___.+$"
    value <- gene <- NULL
    p <- ggplot(df_plt, aes(value, gene)) +
        geom_segment(aes(yend = gene), xend = 0, show.legend = FALSE) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_point(color = "blue") +
        scale_y_discrete(labels = function(x) gsub(reg, "", x)) +
        labs(x = "Loading", y = "Gene")
    if (is_multi) {
        p <- p + facet_wrap(~PC + sample, scales = "free_y", ncol = ncol)
    } else {
        p <- p + facet_wrap(~PC, scales = "free_y", ncol = ncol)
    }
    p
}

.plot_dimdata_bin2d_fun <- function(fun) {
    function(sce, x, y, facet_by = NULL, subset = NULL, bins = 100,
             binwidth = NULL, hex = FALSE, name_true = NULL, name_false = NULL,
             ncol = NULL, ...) {
        args <- list(...)
        if (missing(sce) && "sfe" %in% names(args)) {
            warning("Argument 'sfe' is deprecated. Please use 'sce' instead.")
            sce <- args$sfe
        }
        bin_fun <- if (hex) geom_hex else geom_bin2d
        df <- as.data.frame(fun(sce))
        if (!is.null(facet_by) && !.is_discrete(df[[facet_by]])) {
            warning(facet_by, " is not a categorical variable. Not facetting.")
            facet_by <- NULL
        }
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
        if (!is.null(facet_by)) {
            p <- p + facet_wrap(facet_by, ncol = ncol)
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
#' @param sce A \code{SingleCellExperiment} object.
#' @param x Name of the column in \code{colData} or \code{rowData} to plot on
#'   the x axis of the plot.
#' @param y Name of the column in \code{colData} or \code{rowData} to plot on
#'   the y axis of the plot.
#' @param facet_by Column in \code{colData} or \code{rowData} to facet with.
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
#' @param ncol If facetting, the number of columns of facets, passed to
#'   \code{\link{facet_wrap}}.
#' @importFrom ggplot2 geom_bin2d geom_hex vars
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
    function(sce, feature, fill_by = NULL, facet_by = NULL, subset = NULL,
             bins = 100, binwidth = NULL, scales = "free", ncol = 1,
             position = "stack", ...) {
        args <- list(...)
        if (missing(sce) && "sfe" %in% names(args)) {
            warning("Argument 'sfe' is deprecated. Please use 'sce' instead.")
            sce <- args$sfe
        }
        df <- as.data.frame(fun(sce))[, c(feature, fill_by, facet_by, subset),
                                      drop = FALSE]
        if (!is.null(facet_by) && !.is_discrete(df[[facet_by]])) {
            warning(facet_by, " is not a categorical variable. Not facetting.")
            facet_by <- NULL
        }
        if (!is.null(subset)) df <- df[df[[subset]],]
        p <- ggplot()
        variable <- NULL
        if (length(feature) > 1L) {
            df <- reshape(df, varying = feature, direction = "long",
                          v.names = "values", timevar = "variable",
                          times = feature)
            p <- p +
                geom_histogram(data = df,
                               mapping = aes(!!!syms(c(x = "values", fill = fill_by))),
                               bins = bins,
                               binwidth = binwidth, position = position)
            if (is.null(facet_by)) {
                p <- p + facet_wrap(~ variable, scales = scales, ncol = ncol)
            } else {
                p <- p + facet_grid(rows = vars(variable),
                                    cols = vars(!!!syms(facet_by)),
                                    scales = scales)
            }
        } else {
            p <- p +
                geom_histogram(data = df,
                               mapping = aes(!!!syms(c(x = feature, fill = fill_by))),
                               bins = bins,
                               binwidth = binwidth, position = position)
            if (!is.null(facet_by)) {
                p <- p + facet_wrap(facet_by, ncol = ncol, scales = scales)
            }
        }
        p + scale_fill_manual(values = ditto_colors)
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
#' @param facet_by Column in \code{colData} or \code{rowData} to facet with.
#'   When multiple features are plotted, the features will be in different
#'   facets. In this case, setting \code{facet_by} will call
#'   \code{\link{facet_grid}} so the features are in rows and categories in
#'   \code{facet_by} will be in columns.
#' @param ncol Number of columns in the facetting.
#' @param fill_by Name of a categorical column in \code{colData} or
#'   \code{rowData} to fill the histogram.
#' @param subset Name of a logical column to only plot a subset of the data.
#' @return A ggplot object
#' @seealso plotColDataFreqpoly
#' @importFrom rlang %||% .data
#' @importFrom ggplot2 facet_grid
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
    function(sce, feature, color_by = NULL, subset = NULL, bins = 100,
             binwidth = NULL, linewidth = 1.2,
             scales = "free", ncol = 1, position = "identity") {
        df <- as.data.frame(fun(sce))[, c(feature, color_by, subset), drop = FALSE]
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
