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
    df <- cbind(as.data.frame(rowData(sce)[rownames(loadings), ]),
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
        facet_wrap(~PC, scales = "free_y", ncol = 2) +
        scale_y_discrete(labels = function(x) gsub(reg, "", x)) +
        labs(x = "Loading", y = "Gene")
}
