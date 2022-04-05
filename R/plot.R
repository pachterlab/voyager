# 8. My own plotting function for moran.plot, with ggplot2
# 9. Plotting with divergent palette
# 10. What to do with the image when using geom_sf
# 11. Plot correlograms for multiple genes at once, with error bars (2 sd), as
# in the plot function for spcor, but with ggplot.
# 12. Cluster the correlograms and plot the clusters
# 13. Plot the graphs with spatial coordinates

#' Plot gene expression in space
#'
#' Unlike \code{Seurat} and \code{ggspavis}, plotting functions in this package
#' uses \code{geom_sf} whenever applicable.

plotSFExpression <- function(sfe, colGeometryName, sample_id, color_by,
                             divergent = FALSE,
                             diverge_center = 0, colour_by = color_by,
                             size = 0) {

}


#' Get beginning and end of palette to center a divergent palette
#'
#' The title is self-explanatory.
#'
#' @param values Numeric vector to be colored.
#' @param diverge_center Value to center on, defaults to 0.
#' @return A numeric vector of length 2, the first element is for beginning, and
#' the second for end. The values are between 0 and 1.
#' @export
get_diverge_range <- function(values, diverge_center = 0) {
  rg <- range(values, na.rm = TRUE)
  if (!between(diverge_center, rg[1], rg[2])) {
    stop("diverge_center must be between the minimum and maximum of the metric.")
  }
  rg_centered <- abs(rg - diverge_center)
  if (rg_centered[1] < rg_centered[2]) {
    pal_begin <- (rg_centered[2] - rg_centered[1])/rg_centered[2]/2
    pal_end <- 1
  } else {
    pal_begin <- 0
    pal_end <- 1 - (rg_centered[1] - rg_centered[2])/rg_centered[1]/2
  }
  c(pal_begin, pal_end)
}

#' Use ggplot to plot the moran.plot results
#'
#' Also plots contours showing point density to deal with over-plotting.
#'
#' @param mp Results from spdep::moran.plot.
#' @param var_name Name of the variable to show on the plot. It will be converted
#' to sentence case in the x axis and lower case in the y axis appended after
#' "Spatially lagged".
#' @param cluster string, the column name in mp to plot. The cluster column is
#' added to the mp data frame. Don't use tidyeval. Again, I need a better way to
#' organize the results.
#' @param plot_singletons Logical, whether to plot items that don't have spatial
#' neighbors.
#' @return A ggplot object.
#' @importFrom ggplot2 geom_point aes_string geom_smooth geom_hlive geom_vline
#' geom_density2d scale_shape_manual coord_equal labs
#' @importFrom stringr str_to_sentence str_to_lower
#' @export
moran_ggplot <- function(mp, var_name, cluster = NULL, plot_singletons = TRUE) {
  if (!plot_singletons) {
    mp <- mp[mp$wx > 0,]
  }
  p <- ggplot(mp, aes(x=x, y=wx))
  if (plot_singletons) {
    p <- p +
      geom_point(data = mp[mp$wx == 0,], shape = 21, size = 5, fill = "gray",
                 color = "black")
  }
  if (!is.null(cluster)) {
    pts <- geom_point(aes_string(shape = "is_inf", color = cluster))
  } else {
    pts <- geom_point(aes(shape = is_inf), alpha = 0.7)
  }
  p <- p + pts +
    geom_smooth(formula=y ~ x, method="lm") +
    geom_hline(yintercept=mean(mp$wx), lty=2) +
    geom_vline(xintercept=mean(mp$x), lty=2) +
    geom_density2d() +
    scale_shape_manual(values = c(1, 9)) +
    coord_equal() +
    labs(x = str_to_sentence(var_name),
         y = paste("Spatially lagged", str_to_lower(var_name)),
         shape = "Influential")
  if (!is.null(cluster)) {
    p <- p +
      scale_color_brewer(palette = "Set2")
  }
  p
}

#' Plot uninfluential points from moran.plot as filled contours
#'
#' Just like moran_ggplot, this function plots moran.plot results with ggplot2.
#' However, the uninfluential points are plotted as filled contours and only
#' the influential points are plotted as points
#'
#' @inheritParams moran_ggplot
#' @importFrom ggplot2 geom_density2d_filled scale_fill_viridis_d
#' scale_x_continuous scale_y_continuous expansion
#' @return A ggplot object.
#' @export
moran_ggplot_filled <- function(mp, var_name, plot_singletons = TRUE) {
  if (!plot_singletons) {
    mp <- mp[mp$wx > 0,]
  }
  p <- ggplot(mp, aes(x=x, y=wx))
  if (plot_singletons) {
    p <- p +
      geom_point(data = mp[mp$wx == 0,], shape = 21, size = 5, fill = "gray",
                 color = "black")
  }
  p +
    geom_density2d_filled(show.legend = FALSE) +
    geom_point(data = mp[mp$wx == 0 & mp$is_inf,], shape = 21, size = 5,
               fill = "blue", color = "cornflowerblue") +
    geom_smooth(formula=y ~ x, method="lm") +
    geom_hline(yintercept=mean(mp$wx), lty=2, color = "cornflowerblue") +
    geom_vline(xintercept=mean(mp$x), lty=2, color = "cornflowerblue") +
    geom_point(data=mp[mp$is_inf,], aes(x=x, y=wx), shape=9, color = "cornflowerblue") +
    scale_fill_viridis_d(option = "E") +
    coord_equal() +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion()) +
    labs(x = str_to_sentence(var_name),
         y = paste("Spatially lagged", str_to_lower(var_name)),
         shape = "Influential")
}

#' Plot the elbow plot or scree plot for PCA
#'
#' Apparently, there is no apparent way to plot the PC elbow plot other than
#' extracting the variance explained attribute of the dimred slot, because even
#' the OSCA book makes the elbow plot this way, which I find kind of cumbersome
#' compared to Seurat. So I'm writing this function to make the elbow plot with
#' SCE less cumbersome.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param ndims Number of PCs to plot.
#' @param reduction Name of the dimension reduction to use. It must have an
#' attribute called "percentVar". Defaults to "PCA".
#' @return A ggplot object. The y axis is percentage of variance explained.
#' @export
ElbowPlot <- function(sce, ndims = 20, reduction = "PCA") {
  percent_var <- attr(reducedDim(sce, reduction), "percentVar")
  inds <- seq_len(ndims)
  df <- data.frame(PC = inds,
                   pct_var = percent_var[inds])
  ggplot(df, aes(PC, pct_var)) +
    geom_point() +
    labs(x = "PC", y = "Variance explained (%)")
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
#' @param show_symbol Logical; if the row names of the matrix are Ensembl accessions,
#' indicate whether to show more human readable gene symbols in the plot instead.
#' Ignored if the column specified in \code{symbol_col} is absent from rowData.
#' @param symbol_col If the row names of the gene expression matrix are Ensembl accessions
#' to avoid ambiguity in analysis, specify
#' @param balanced Return an equal number of genes with + and - scores. If
#'   FALSE, returns the top genes ranked by the scores absolute values.
#' @param ncol Number of columns in the facetted plot.
#' @return A ggplot object. Loadings for different PCs are plotted in different
#' facets so one ggplot object is returned.

plotDimLoadings <- function(sce, dims = 1:4, nfeatures = 10,
                            show_symbol = TRUE, symbol_col = "symbol",
                            reduction = "PCA",
                            balanced = TRUE, ncol = 2) {
  loadings <- attr(reducedDim(sce, reduction), "rotation")
  is_ensembl <- all(str_detect(rownames(sce), "^ENS"))
  if (!is_ensembl) show_symbol <- FALSE # I mean, irrelevant
  loading_cols <- paste0("PC", dims)
  df <- cbind(as.data.frame(rowData(sce)[rownames(loadings),]), loadings[,loading_cols])
  df <- df %>%
    pivot_longer(starts_with("PC"), names_to = "PC") %>%
    group_by(PC)
  if (!symbol_col %in% names(df) || !show_symbol) {
    df$gene_show <- rownames(sce)
  } else {
    df$gene_show <- df[[symbol_col]]
  }
  if (balanced) {
    n2 <- floor(nfeatures/2)
    df_plt_plus <- df %>%
      slice_max(value, n = n2)
    df_plt_minus <- df %>%
      slice_max(-value, n = n2)
    df_plt <- rbind(df_plt_plus, df_plt_minus)
  } else {
    df_plt <- df %>%
      slice_max(abs(value), n = nfeatures)
  }
  df_plt <- df_plt %>%
    ungroup() %>%
    mutate(PC = fct_relevel(PC, loading_cols),
           gene = tidytext::reorder_within(gene_show, value, PC))
  ggplot(df_plt, aes(value, gene)) +
    geom_segment(aes(yend = gene), xend = 0, show.legend = FALSE) +
    geom_point(color = "blue") +
    geom_vline(xintercept = 0, linetype = 2) +
    facet_wrap(~ PC, scales = "free_y", ncol = 2) +
    tidytext::scale_y_reordered()
}
