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

.get_applicable <- function(type, fixed) {
  fixed <- .fill_defaults(fixed)
  if (type %in% c("POINT", "MULTIPOINT")) {
    names_use <- c("size", "shape", "alpha", "color")
    if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 1
    shape <- fixed[["shape"]]
    if (!is.null(shape) && shape > 20L) {
      names_use <- c(names_use, "fill")
    }
  } else if (type %in% c("LINESTRING", "MULTILINESTRING")) {
    if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 1
    names_use <- c("size", "linetype", "alpha", "color")
  } else {
    # i.e. polygons
    names_use <- c("size", "linetype", "fill", "color", "alpha")
  }
  fixed_applicable <- .drop_null_list(fixed[names_use])
  fixed_applicable
}

.get_pal <- function(df, feature_aes, option, divergent, diverge_center) {
  feature_aes <- feature_aes[names(feature_aes) %in% c("fill", "color")]
  if (length(feature_aes)) {
    # color_by is set to NULL if fill_by is applicable and present
    m <- df[[unlist(feature_aes)]]
    .aes <- names(feature_aes)
  } else {
    return(NULL)
  }
  if (.is_discrete(m)) {
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
      pal_fun <- switch(.aes, fill = scale_fill_scico,
                        color = scale_color_scico)
      pal <- pal_fun(palette = .pal, begin = pal_begin,
                     end = pal_end, na.value = "gray")
    } else {
      pal_fun <- switch(.aes, fill = scale_fill_distiller,
                        color = scale_color_distiller)
      .pal <- switch(option, "Blues", "PuRd")
      pal <- pal_fun(na.value = "gray", palette = .pal, direction = 1)
    }
  }
  pal
}

.fill_defaults <- function(fixed) {
  defaults <- list(size = 0, shape = 16,
                   linetype = 1, alpha = 1, color = "black", fill = "gray70",
                   divergent = FALSE, diverge_center = NULL)
  fill <- defaults[setdiff(names(fixed), names(defaults))]
  .drop_null_list(c(fixed, fill))
}

#' @importFrom sf st_is st_drop_geometry st_geometry_type
#' @importFrom ggplot2 ggplot aes_string geom_sf scale_fill_manual
#' scale_color_manual scale_fill_distiller scale_color_distiller
#' @importFrom scico scale_fill_scico scale_color_scico
#' @importFrom ggnewscale new_scale_color
.plot_var_sf <- function(df, annot_df, type, type_annot, feature_aes, feature_fixed,
                         annot_aes, annot_fixed, divergent, diverge_center) {
  # Add annotGeometry if present
  if (!is.null(annot_df)) {
    annot_fixed <- .get_applicable(type_annot, annot_fixed)
    annot_fixed <- annot_fixed[setdiff(names(annot_fixed), names(annot_aes))]
    aes_annot <- do.call(aes_string, annot_aes)
    geom_annot <- do.call(geom_sf, c(list(mapping = aes_annot, data = annot_df),
                                     annot_fixed))
    pal_annot <- .get_pal(annot_df, annot_aes, 2, annot_fixed$divergent,
                          annot_fixed$diverge_center)
  }

  p <- ggplot()
  # Filled polygon annotations go beneath feature plot
  is_annot_filled <- !is.null(annot_df) &&
    ("fill" %in% names(c(annot_aes, annot_fixed))) &&
    (st_is(annot_df, "POLYGON") || st_is(annot_df, "MULTIPOLYGON"))
  if ("fill" %in% names(annot_fixed))
    is_annot_filled <- is_annot_filled && !is.na(annot_fixed[["fill"]])
  if (is_annot_filled) {
    p <- p + geom_annot
    if (!is.null(pal_annot)) p <- p + pal_annot
  }

  feature_fixed <- .get_applicable(type, feature_fixed)
  feature_fixed <- feature_fixed[setdiff(names(feature_fixed), names(feature_aes))]

  if ("fill" %in% names(feature_aes) && is_annot_filled)
    p <- p + new_scale_fill()
  aes_use <- do.call(aes_string, feature_aes)
  geom_use <- do.call(geom_sf, c(list(mapping = aes_use, data = df),
                                 feature_fixed))
  p <- p + geom_use

  # Palette
  pal <- .get_pal(df, feature_aes, 1, divergent, diverge_center)
  if (!is.null(pal)) p <- p + pal

  # Line and point annotations go above feature plot
  if (!is.null(annot_df) && !is_annot_filled) {
    if (!is.null(pal_annot) && "color" %in% names(feature_aes)) {
      p <- p + new_scale_color()
    }
    p <- p + geom_annot
    if (!is.null(pal_annot)) {
      p <- p + pal_annot
    }
  }
  p
}

.get_feature_aes <- function(m, type, aes_spec, shape) {
  if (type %in% c("POINT", "MULTIPOINT")) {
    if (aes_spec == "fill" && shape < 20) aes_spec <- "color"
    if (aes_spec == "linetype") aes_spec <- "shape"
  }
  if (type %in% c("LINESTRING", "MULTILINESTRING")) {
    if (aes_spec == "fill") aes_spec <- "color"
    if (aes_spec == "shape") aes_spec <- "linetype"
  }
  if (!.is_discrete(m) && aes_spec %in% c("shape", "linetype")) {
    stop("Shape and linetype are only applicable to discrete variables.")
  }
  aes_spec
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
#' @param divergent Logical, whether a divergent palette should be used.
#' @param diverge_center If \code{divergent = TRUE}, the center from which the
#'   palette should diverge. If \code{NULL}, then not centering.
#' @param size Fixed size of points or width of lines, including outlines of
#'   polygons. For polygons, this defaults to 0, meaning no outlines. For points
#'   and lines, this defaults to 1. Ignored if \code{size_by} is specified.
#' @param shape Fixed shape of points, ignored if \code{shape_by} is specified
#'   and applicable.
#' @param linetype Fixed line type, ignored if \code{linetype_by} is specified
#'   and applicable.
#' @param color Fixed color for \code{colGeometry} if \code{color_by} is not
#'   specified or not applicable, or for \code{annotGeometry} if
#'   \code{annot_color_by} is not specified or not applicable.
#' @param fill Similar to \code{color}, but for fill.
#' @param alpha Transparency.
#' @param annotGeometryName Name of a \code{annotGeometry} of the SFE object, to
#'   annotate the gene expression plot.
#' @param annot_aes A named list of plotting parameters for the annotation sf
#'   data frame. The names are which geom (as in ggplot2, such as color and
#'   fill), and the values are column names in the annotation sf data frame.
#'   Tidyeval is NOT supported.
#' @param annot_fixed Similar to \code{annot_aes}, but for fixed aesthetic
#'   settings, such as \code{color = "gray"}. The defaults are the same as the
#'   relevant defaults for this function.
#' @param ncol Number of columns if plotting multiple features. Defaults to
#'   \code{NULL}, which means using the same logic as \code{facet_wrap}, which
#'   is used by \code{patchwork}'s \code{\link{wrap_plots}} by default.
#' @param aes_use Aesthetic to use for discrete variables. For continuous
#'   variables, it's always "fill" for polygons and point shapes 21-25. For
#'   discrete variables, it can be fill, color, shape, or linetype, whenever
#'   applicable. The specified value will be changed to the applicable equivalent.
#'   For example, if the geometry is point but "linetype" is specified, then
#'   "shaped" will be used instead.
#' @param ... Other arguments passed to \code{\link{wrap_plots}}.
#' @importFrom patchwork wrap_plots
#' @importFrom stats setNames
#' @importMethodsFrom Matrix t
plotSpatialFeature <- function(sfe, colGeometryName, features, sample_id = NULL,
                               ncol = NULL, annotGeometryName = NULL,
                               annot_aes = list(), annot_fixed = list(),
                               exprs_values = "logcounts",
                               aes_use = c("fill", "color", "shape", "linetype"),
                               divergent = FALSE,
                               diverge_center = NULL, only_plot_expressed = FALSE,
                               size = 0, shape = 16, linetype = 1, alpha = 1,
                               color = "black", fill = "gray70", ...) {
  aes_use <- match.arg(aes_use)
  features_list <- .check_features(sfe, features, colGeometryName)
  values <- list()
  if (is.null(sample_id))
    sample_id_ind <- rep(TRUE, ncol(sfe))
  else
    sample_id_ind <- colData(sfe)$sample_id %in% sample_id
  if (!is.null(features_list[["assay"]])) {
    values_assay <- assay(sfe, exprs_values)[features_list[["assay"]],
                                             sample_id_ind, drop = FALSE]
    values_assay <- as.data.frame(as.matrix(t(values_assay)))
    values[["assay"]] <- values_assay
  }
  if (!is.null(features_list[["coldata"]]))
    values[["coldata"]] <- as.data.frame(colData(sfe)[sample_id_ind,
                                                      features_list[["coldata"]],
                                                      drop = FALSE])
  if (length(values) > 1L) values <- cbind(values$assay, values$coldata)

  df <- colGeometry(sfe, colGeometryName, sample_id = sample_id)
  # Will use separate ggplots for each feature so each can have its own color scale
  if (!is.null(annotGeometryName)) {
    annot_df <- annotGeometry(sfe, annotGeometryName, sample_id)
    type_annot <- st_geometry_type(annot_df, by_geometry = FALSE)
  }
  else annot_df <- NULL; type_annot <- NULL
  feature_fixed <- list(size = size, shape = shape, linetype = linetype,
                        alpha = alpha, color = color, fill = fill)
  type <- st_geometry_type(df, by_geometry = FALSE)
  plots <- lapply(names(values), function(n) {
    df[[n]] <- values[[n]]
    feature_aes_name <- .get_feature_aes(df[[n]], type, aes_use, shape)
    feature_aes <- setNames(list(n), feature_aes_name)
    .plot_var_sf(df, annot_df, type, type_annot, feature_aes, feature_fixed,
                 annot_aes, annot_fixed, divergent, diverge_center)
  })
  if (length(plots) > 1L) {
    out <- wrap_plots(plots, ncol = ncol, ...)
  } else {
    out <- plots[[1]]
  }
  return(out)
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
#' @importFrom ggplot2 geom_point aes_string geom_smooth geom_hline geom_vline
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
