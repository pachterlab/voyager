# 8. My own plotting function for moran.plot, with ggplot2
# 9. Plotting with divergent palette
# 10. What to do with the image when using geom_sf
# 11. Plot correlograms for multiple genes at once, with error bars (2 sd), as
# in the plot function for spcor, but with ggplot.
# 12. Cluster the correlograms and plot the clusters
# 13. Plot the graphs with spatial coordinates

#' Get beginning and end of palette to center a divergent palette
#'
#' The title is self-explanatory.
#'
#' @param values Numeric vector to be colored.
#' @param diverge_center Value to center on, defaults to 0.
#' @return A numeric vector of length 2, the first element is for beginning, and
#' the second for end. The values are between 0 and 1.
#' @export
getDivergeRange <- function(values, diverge_center = 0) {
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

.get_applicable <- function(df, fill_by, color_by, shape_by, linetype_by,
                            size_by, size, shape, linetype, alpha, color,
                            fill) {
  if (st_is(df, "POINT") || st_is(df, "MULTIPOINT")) {
    aes_applicable <- list(geometry = "geometry",
                           color = color_by, shape = shape_by,
                           size = size_by)
    if (isTRUE(all.equal(0, size))) size <- 1
    fixed_applicable <- list(size = size, shape = shape, alpha = alpha,
                             color = color)
    if (!is.null(shape) && shape > 20L) {
      if (!is.null(fill_by)) aes_applicable[["color"]] <- NULL
      aes_applicable <- c(aes_applicable, list(fill = fill_by))
      fixed_applicable <- c(fixed_applicable, list(fill = fill))
    }
  } else if (st_is(df, "LINESTRING") || st_is(df, "MULTILINESTRING")) {
    aes_applicable <- list(geometry = "geometry", linetype = linetype_by,
                           color = color_by, shape = shape_by, size = size_by)
    if (isTRUE(all.equal(0, size))) size <- 1
    fixed_applicable <- list(size = size, linetype = linetype, alpha = alpha,
                             color = color)
  } else {
    if (!is.null(fill_by)) color_by <- NULL
    aes_applicable <- list(geometry = "geometry", fill = fill_by,
                           color = color_by, size = size_by,
                           linetype = linetype_by)
    fixed_applicable <- list(size = size, linetype = linetype, alpha = alpha,
                             color = color, fill = fill)
  }
  aes_applicable <- .drop_null_list(aes_applicable)
  fixed_applicable <- .drop_null_list(fixed_applicable)
  fixed_applicable <- fixed_applicable[setdiff(names(fixed_applicable),
                                               names(aes_applicable))]
  list(aes_applicable, fixed_applicable)
}

.get_pal <- function(df, applicable, option, divergent, diverge_center) {
  cols_check_names <- unlist(applicable)
  cols_check_names <- cols_check_names[names(cols_check_names) %in% c("fill", "color")]
  if (length(cols_check_names)) {
    # color_by is set to NULL if fill_by is applicable and present
    m <- st_drop_geometry(df)[,cols_check_names]
    is_discrete <- function(m) is.character(m) | is.factor(m) | is.logical(m)
    .aes <- names(cols_check_names)
  } else {
    return(NULL)
  }
  if (is_discrete) {
    .pal <- switch (option, Voyager::ditto_colors, rev(Voyager::ditto_colors))
    pal_fun <- switch (.aes, fill = scale_fill_manual, color = scale_color_manual)
    pal <- pal_fun(values = .pal, na.value = "gray")
  } else {
    if (divergent) {
      if (!is.null(diverge_center)) {
        r <- df[[aes_applicable[["fill"]]]]
        pal_range <- getDivergeRange(r, diverge_center)
        pal_begin <- pal_range[1]
        pal_end <- pal_range[2]
      } else {
        pal_begin <- 0
        pal_end <- 1
      }
      .pal <- switch(option, "roma", "bam")
      pal_fun <- switch (.aes, fill = scale_fill_scico,
                         color = scale_color_scico)
      pal <- pal_fun(palette = .pal, begin = pal_begin,
                     end = pal_end, na.value = "gray")
    } else {
      pal_fun <- switch(.aes, fill = scale_fill_brewer,
                        color = scale_color_brewer)
      .pal <- switch(option, "Blues", "PuRd")
      pal <- pal_fun(na.value = "gray", palette = .pal, direction = 1)
    }
  }
  pal
}

.annot_defaults <- function(annot_params) {
  defaults <- list(fill_by = NULL, color_by = NULL, shape_by = NULL,
                   linetype_by = NULL, size_by = NULL, size = 0, shape == 16,
                   linetype = 1, alpha = 1, color = "black", fill = "gray70",
                   divergent = FALSE, diverge_center = NULL)
  if (!is.null(annot_params$fill_by)) annot_params$color_by <- NULL
  fill <- defaults[setdiff(names(annot_params), names(defaults))]
  .drop_null_list(c(annot_params, fill))
}

#' @importFrom sf st_is st_drop_geometry
#' @importFrom ggplot2 ggplot aes_string geom_sf scale_fill_manual
#' scale_color_manual scale_fill_brewer scale_color_brewer
#' @importFrom scico scale_fill_scico scale_color_scico
#' @importFrom ggnewscale new_scale_color
.plot_var_sf <- function(df, fill_by, color_by, shape_by, linetype_by, size_by,
                         annot_df, annot_params, divergent, diverge_center,
                         only_plot_expressed, size, shape, linetype, alpha,
                         color, fill) {
  # Add annotGeometry if present
  if (!is.null(annot_df)) {
    annot_params <- .annot_defaults(annot_params)
    .by <- grepl("_by$", names(annot_aes))
    annot_fixed <- annot_params[!.by]
    annot_aes <- annot_params[.by]
    aes_annot <- do.call(aes_string, annot_aes)
    geom_annot <- do.call(geom_sf, c(list(mapping = aes_use, data = annot_df),
                                     annot_fixed))
    pal_annot <- .get_pal(annot_df, annot_aes, 2, annot_params$divergent,
                          annot_params$diverge_center)
  }

  p <- ggplot()
  # Polygon annotations go beneath feature plot
  is_annot_polygon <- !is.null(annot_df) && (st_is(annot_df, "POLYGON") || st_is(annot_df, "MULTIPOLYGON"))
  if (is_annot_polygon) {
    p <- p + geom_annot
    if (!is.null(pal_annot)) p <- p + pal_annot
  }

  applicable <- .get_applicable(df, fill_by, color_by, shape_by, linetype_by,
                                size_by, size, shape, linetype, alpha, color,
                                fill)
  if (only_plot_expressed) {
    col_filter <- unlist(applicable[[1]][names(applicable[[1]]) %in% c("fill", "color")])
    if (all(df[[col_filter]] >= 0)) {
      df <- df[df[[col_filter]] > 0,]
    }
  }

  if ("fill" %in% names(applicable[[1]]) && is_annot_polygon)
    p <- p + new_scale_fill()
  aes_use <- do.call(aes_string, applicable[[1]])
  geom_use <- do.call(geom_sf, c(list(mapping = aes_use, data = df),
                                 applicable[[2]]))
  p <- p + geom_use

  # Palette
  pal <- .get_pal(df, applicable[[1]], 1, divergent, diverge_center)
  if (!is.null(pal)) p <- p + pal

  # Line and point annotations go above feature plot
  if (!is.null(annot_df) && !is_annot_polygon) {
    if (!is.null(pal_annot) && "color" %in% names(applicable[[1]])) {
      p <- p + new_scale_color()
    }
    p <- p + geom_use
    if (!is.null(pal_annot)) {
      p <- p + pal_annot
    }
  }
  p
}

#' Plot gene expression in space
#'
#' Unlike \code{Seurat} and \code{ggspavis}, plotting functions in this package
#' uses \code{geom_sf} whenever applicable.
#'
#' In the documentation of this function, a "feature" can be a gene (or whatever
#' entity that corresponds to rows of the gene count matrix), a column in
#' \code{colData}, or a column in the \code{colGeometry} \code{sf} data frame
#' specified in the \code{colGeometryName} argument.
#'
#' For continuous variables, the Blues palette from colorbrewer is used if
#' \code{divergent = FALSE}, and the roma palette from the \code{scico} package
#' if \code{divergent = TRUE}. For discrete variables, the \code{dittoSeq}
#' palette is used. The defaults are colorblind friendly. For annotation, the
#' PuRd colorbrewer palette is used for continuous variables and the other end
#' of the \code{dittoSeq} palette is used for discrete variables.
#'
#' @inheritParams calculateMoransI
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param fill_by Feature to fill the polygons or point shapes that supports
#'   fill. The polygons will not be filled if \code{NULL}.
#' @param color_by Feature to color the points or lines, including outlines of
#'   polygons. For polygons, this is ignored if \code{fill_by} is specified to
#'   avoid an overly garish and hard to read plot.
#' @param shape_by Feature for shape of points, only applicable if the
#'   \code{colGeometry} from \code{colGeometryName} is of type POINT or
#'   MULTIPOINT.
#' @param linetype_by Feature for line type, only applicable for LINESTRING,
#'   MULTILINESTRING, POLYGON, and MULTIPOLYGON.
#' @param divergent Logical, whether a divergent palette should be used.
#' @param diverge_center If \code{divergent = TRUE}, the center from which the
#'   palette should diverge. If \code{NULL}, then not centering.
#' @param only_plot_expressed Logical, if \code{TRUE}, and if all values are
#'   non-negative, then geometries with value 0 are not plotted.
#' @param colour_by Same as color_by.
#' @param size Fixed size of points or width of lines, including outlines of
#'   polygons. For polygons, this defaults to 0, meaning no outlines. For points
#'   and lines, this defaults to 1. Ignored if \code{size_by} is specified.
#' @param shape Fixed shape of points, ignored if \code{shape_by} is specified
#'   and applicable.
#' @param linetype Fixed line type, ignored if \code{linetype_by} is specified
#'   and applicable.
#' @param color Fixed color for \code{colGeometry} if \code{color_by} is not
#' specified or not applicable, or for \code{annotGeometry} if \code{annot_color_by}
#' is not specified or not applicable.
#' @param fill Similar to \code{color}, but for fill.
#' @param alpha Transparency.
#' @param annotGeometryName Name of a \code{annotGeometry} of the SFE object, to
#'   annotate the gene expression plot.
#' @param annot_color_by Same as \code{color_by}, but for the annotation when
#'   \code{annotGeometryName} is specified.
#' @param annot_colour_by Same as \code{annot_color_by}.
plotSpatialFeature <- function(sfe, colGeometryName, features, sample_id,
                               fill_by = NULL, color_by = NULL, shape_by = NULL,
                               linetype_by = NULL, size_by = NULL,
                               annotGeometryName = NULL, annot_aes = list(),
                               exprs_values = "logcounts", divergent = FALSE,
                               diverge_center = NULL, only_plot_expressed = FALSE,
                               colour_by = color_by, size = 0,
                               shape = 16, linetype = 1, alpha = 1, color = "black",
                               fill = "gray70") {
  features_list <- .check_features(sfe, features, colGeometryName)
  features_use <- assay(sfe, exprs_values)[features, colData(sfe)$sample_id %in% sample_id]

  df <- colGeometry(sfe, colGeometryName, sample_id = sample_id)
  if (length(features) == 1L) {
    df[[features]] <- features_use
  } else {
    features_use <- t(as.matrix(features_use))
    features_use <- as.data.frame(features_use)
    df <- cbind(df, features_use)
  }
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
