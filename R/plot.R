# 10. What to do with the image when using geom_sf
# 11. Plot correlograms for multiple genes at once, with error bars (2 sd), as
# in the plot function for spcor, but with ggplot.
# 12. Cluster the correlograms and plot the clusters
# 14. Plot MoranMC results with ggplot2
# Shall I also allow users to plot dimension reductions as features?
# For example, plotting PC1 in space, as opposed to MULTISPATI PC1. I think I'll
# do that, not only for plotting functions, but also for the metrics.
# To do:
# 1. I kind of find it annoying to type in colGeometryName and colGraphName.
# Make those optional when there's only one anyway. Generalize .check_sample_id
# 2. reverse_y? It's kind of hard to do that with sf.

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
  if (!(diverge_center >= rg[1] && diverge_center <= rg[2])) {
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
    if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 0.5
    shape <- fixed[["shape"]]
    if (!is.null(shape) && shape > 20L) {
      names_use <- c(names_use, "fill")
    }
  } else if (type %in% c("LINESTRING", "MULTILINESTRING")) {
    if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 0.5
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
        r <- df[[feature_aes[[.aes]]]]
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
                   linetype = 1, alpha = 1, color = "black", fill = "gray80",
                   divergent = FALSE, diverge_center = NULL)
  fill <- defaults[setdiff(names(defaults), names(fixed))]
  out <- .drop_null_list(c(fixed, fill))
  if (is.na(out$fill)) out$size <- 0.5
  out
}

#' @importFrom sf st_drop_geometry st_geometry_type
#' @importFrom ggplot2 ggplot aes_string geom_sf scale_fill_manual
#' scale_color_manual scale_fill_distiller scale_color_distiller geom_polygon
#' geom_segment stat_density2d
#' @importFrom scico scale_fill_scico scale_color_scico
#' @importFrom ggnewscale new_scale_color new_scale_fill
.plot_var_sf <- function(df, annot_df, type, type_annot, feature_aes, feature_fixed,
                         annot_aes, annot_fixed, divergent, diverge_center,
                         annot_divergent, annot_diverge_center) {
  # Add annotGeometry if present
  if (!is.null(annot_df)) {
    annot_fixed <- .get_applicable(type_annot, annot_fixed)
    annot_fixed <- annot_fixed[setdiff(names(annot_fixed), names(annot_aes))]
    if ("color" %in% names(annot_aes) && annot_fixed$size == 0)
      annot_fixed$size <- 0.5
    aes_annot <- do.call(aes_string, annot_aes)
    geom_annot <- do.call(geom_sf, c(list(mapping = aes_annot, data = annot_df),
                                     annot_fixed))
    pal_annot <- .get_pal(annot_df, annot_aes, 2, annot_divergent, annot_diverge_center)
  }

  p <- ggplot()
  # Filled polygon annotations go beneath feature plot
  is_annot_filled <- !is.null(annot_df) &&
    ("fill" %in% names(c(annot_aes, annot_fixed))) &&
    type_annot %in% c("POLYGON", "MULTIPOLYGON")
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
  if (type %in% c("POLYGON", "MULTIPOLYGON") || shape >= 20) {
    aes_spec <- "fill"
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
#' @param features Features to plot, must be in rownames of the gene count
#'   matrix, colnames of colData or a colGeometry.
#' @param divergent Logical, whether a divergent palette should be used.
#' @param diverge_center If \code{divergent = TRUE}, the center from which the
#'   palette should diverge. If \code{NULL}, then not centering.
#' @param size Fixed size of points or width of lines, including outlines of
#'   polygons. For polygons, this defaults to 0, meaning no outlines. For points
#'   and lines, this defaults to 0.5. Ignored if \code{size_by} is specified.
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
#'   applicable. The specified value will be changed to the applicable
#'   equivalent. For example, if the geometry is point but "linetype" is
#'   specified, then "shaped" will be used instead.
#' @param annot_divergent Just as \code{divergent}, but for the annotGeometry in
#'   case it's different.
#' @param annot_diverge_center Just as \code{diverge_center}, but for the
#'   annotGeometry in case it's different.
#' @param only_plot_expressed Logical, whether to only plot values > 0. This
#' argument is only used when all values are non-negative.
#' @param ... Other arguments passed to \code{\link{wrap_plots}}.
#' @importFrom patchwork wrap_plots
#' @importFrom stats setNames
#' @importMethodsFrom Matrix t
plotSpatialFeature <- function(sfe, colGeometryName, features, sample_id = NULL,
                               ncol = NULL, annotGeometryName = NULL,
                               annot_aes = list(), annot_fixed = list(),
                               exprs_values = "logcounts",
                               aes_use = c("fill", "color", "shape", "linetype"),
                               divergent = FALSE, diverge_center = NULL,
                               annot_divergent = FALSE,
                               annot_diverge_center = NULL,
                               only_plot_expressed = FALSE,
                               size = 0, shape = 16, linetype = 1, alpha = 1,
                               color = "black", fill = "gray80", ...) {
  aes_use <- match.arg(aes_use)
  sample_id <- .check_sample_id(sfe, sample_id)
  values <- .get_feature_values(sfe, features, sample_id,
                                colGeometryName = colGeometryName,
                                exprs_values = exprs_values)
  df <- colGeometry(sfe, colGeometryName, sample_id = sample_id)
  # Will use separate ggplots for each feature so each can have its own color scale
  if (!is.null(annotGeometryName)) {
    annot_df <- annotGeometry(sfe, annotGeometryName, sample_id)
    type_annot <- st_geometry_type(annot_df, by_geometry = FALSE)
  }
  else {
    annot_df <- NULL
    type_annot <- NULL
  }
  feature_fixed <- list(size = size, shape = shape, linetype = linetype,
                        alpha = alpha, color = color, fill = fill)
  type <- st_geometry_type(df, by_geometry = FALSE)
  plots <- lapply(names(values), function(n) {
    df[[n]] <- values[[n]]
    feature_aes_name <- .get_feature_aes(df[[n]], type, aes_use, shape)
    feature_aes <- setNames(list(n), feature_aes_name)
    .plot_var_sf(df, annot_df, type, type_annot, feature_aes, feature_fixed,
                 annot_aes, annot_fixed, divergent, diverge_center,
                 annot_divergent, annot_diverge_center)
  })
  if (length(plots) > 1L) {
    out <- wrap_plots(plots, ncol = ncol, ...)
  } else {
    out <- plots[[1]]
  }
  return(out)
}

#' @importFrom ggplot2 geom_line
.plot_graph <- function(listw, coords, geometry = NULL, segment_size = 0.5,
                        geometry_size = 0.5) {
  cardnb <- card(listw$neighbours)
  n <- length(listw$neighbours)
  df <- data.frame(i = rep(1:n, cardnb),
                   j = unlist(listw$neighbours))
  df$x <- coords[,1][df$i]
  df$y <- coords[,2][df$i]
  df$x_end <- coords[,1][df$j]
  df$y_end <- coords[,2][df$j]
  p <- ggplot(df)
  X <- Y <- group <- x <- y <- x_end <- y_end <- NULL
  if (!is.null(geometry)) {
    # Burning question: Shall I use geom_sf or just get the coordinates?
    # Maybe for now I'll just get the coordinates.
    coord_geom <- as.data.frame(st_coordinates(geometry))
    names(coord_geom)[ncol(coord_geom)] <- "group"
    geom_use <- switch (as.character(st_geometry_type(geometry, by_geometry = FALSE)),
                        POINT = geom_point(data = coord_geom, aes(X, Y), size = geometry_size,
                                           color = "gray70"),
                        LINESTRING = geom_line(data = coord_geom, aes(X, Y, group = group),
                                               size = geometry_size, color = "gray70"),
                        POLYGON = geom_polygon(data = coord_geom, aes(X, Y, group = group),
                                               size = geometry_size, color = "gray70", fill = NA)
    )
    p <- p + geom_use
  } else {
    p <- p +
      geom_point(aes(x, y), size = geometry_size)
  }
  p +
    geom_segment(aes(x, y, xend = x_end, yend = y_end), size = segment_size) +
    coord_equal() +
    labs(x = NULL, y = NULL)
}

#' Plot spatial graphs
#'
#' A ggplot version of \code{spdep::plot.nb}, reducing boilerplate for SFE
#' objects.
#'
#' @inheritParams plotSpatialFeature
#' @param colGraphName Name of graph associated with columns of the gene count
#'   matrix to be plotted.
#' @param segment_size Thickness of the segments that represent graph edges.
#' @param geometry_size Point size (for POINT geometries) or line thickness (for
#'   LINESTRING and POLYGON) to plot the geometry in the background.
#' @param annotGraphName Name of the annotation graph to plot.
#' @param annotGeometryName Name of the \code{annotGeometry}, which is
#'   associated with the graph specified with \code{annotGraphName}, for spatial
#'   coordinates of the graph nodes and for context.
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom spdep card
#' @importFrom sf st_coordinates st_centroid st_geometry
#' @return A ggplot2 object.
#' @export
plotColGraph <- function(sfe, colGraphName, colGeometryName = NULL, sample_id = NULL,
                         segment_size = 0.5, geometry_size = 0.5) {
  sample_id <- .check_sample_id(sfe, sample_id)
  g <- colGraph(sfe, colGraphName, sample_id)
  if (is.null(g)) stop("Graph ", colGraphName, " not found.")
  coords <- spatialCoords(sfe)[colData(sfe)$sample_id == sample_id,]
  geometry <- if (is.null(colGeometryName)) NULL else
    colGeometry(sfe, colGeometryName, sample_id)
  .plot_graph(g, coords, geometry, segment_size, geometry_size)
}

#' @rdname plotColGraph
#' @export
plotAnnotGraph <- function(sfe, annotGraphName, annotGeometryName,
                           sample_id = NULL, segment_size = 0.5,
                           geometry_size = 0.5) {
  sample_id <- .check_sample_id(sfe, sample_id)
  g <- annotGraph(sfe, annotGraphName, sample_id)
  ag <- annotGeometry(sfe, annotGeometryName, sample_id)
  ag_type <- st_geometry_type(ag, by_geometry = FALSE)
  if (ag_type == "POINT") {
    coords <- st_coordinates(ag)
  } else {
    coords <- st_coordinates(st_centroid(st_geometry(ag)))
  }
  .plot_graph(g, coords, ag, segment_size, geometry_size)
}

.moran_ggplot <- function(mp, feature, is_singleton, contour_color = "cyan",
                          color_by = NULL, plot_singletons = TRUE,
                          divergent = FALSE, diverge_center = NULL, ...) {
  if (!plot_singletons) {
    mp <- mp[!is_singleton,]
  }
  x <- wx <- is_inf <- NULL
  if (all(!is_singleton) && plot_singletons) plot_singletons <- FALSE
  p <- ggplot(mp, aes(x=x, y=wx))
  if (plot_singletons) {
    # Need to use listw to check for singletons.
    p <- p +
      geom_point(data = mp[is_singleton,], shape = 21, size = 5, fill = "gray",
                 color = "black")
  }
  if (!is.null(color_by)) {
    pal <- .get_pal(mp, list(color = color_by), option = 1,
                    divergent = divergent, diverge_center = diverge_center)
    p <- p + pal
    pts <- geom_point(aes_string(shape = "is_inf", color = color_by))
  } else {
    pts <- geom_point(aes(shape = is_inf), alpha = 0.5)
  }
  p <- p + pts
  # stat_density2d doesn't work when there're too few points
  # Unlikely in real data, but just in case
  # The error doesn't show up until the plot is built.
  p_test <- tryCatch(ggplot_build(p + stat_density2d(...)),
                     error = function(e) {
                       warning("Too few points for stat_density2d, not plotting contours.")
                     },
                     warning = function(w) {
                       warning("Too few points for stat_density2d, not plotting contours.")
                     })
  if (is(p_test, "ggplot_built"))
    p <- p + geom_density2d(color = contour_color, ...)
  p <- p +
    geom_smooth(formula=y ~ x, method="lm") +
    geom_hline(yintercept=mean(mp$wx), lty=2, color = "gray") +
    geom_vline(xintercept=mean(mp$x), lty=2, color = "gray") +
    scale_shape_manual(values = c(1, 9)) +
    coord_equal() +
    labs(x = feature,
         y = paste("Spatially lagged", feature),
         shape = "Influential")
  p
}

.moran_ggplot_filled <- function(mp, feature, is_singleton, color_by = NULL,
                                 plot_singletons = TRUE, divergent = FALSE,
                                 diverge_center = NULL, ...) {
  if (!plot_singletons) {
    mp <- mp[!is_singleton,]
  }
  x <- wx <- is_inf <- NULL
  p <- ggplot(mp, aes(x=x, y=wx)) +
    geom_density2d_filled(show.legend = FALSE, ...)
  if (plot_singletons) {
    p <- p +
      geom_point(data = mp[is_singleton & mp$is_inf,], shape = 21, size = 5,
                 fill = "blue", color = "cornflowerblue")
  }
  mp_inf <- mp[mp$is_inf,]
  if (!is.null(color_by)) {
    pal <- .get_pal(mp_inf, list(color = color_by), option = 1,
                    divergent = divergent, diverge_center = diverge_center)
    p <- p + pal
    pts <- geom_point(data = mp_inf, aes_string(color = color_by),
                      shape = 9)
  } else {
    pts <- geom_point(data = mp_inf, shape = 9, color = "cornflowerblue")
  }
  p +
    geom_smooth(formula=y ~ x, method="lm") +
    geom_hline(yintercept=mean(mp$wx), lty=2, color = "cornflowerblue") +
    geom_vline(xintercept=mean(mp$x), lty=2, color = "cornflowerblue") +
    pts +
    scale_fill_viridis_d(option = "E") +
    coord_equal() +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion()) +
    labs(x = feature,
         y = paste("Spatially lagged", feature),
         shape = "Influential")
}

#' Use ggplot to plot the moran.plot results
#'
#' This function uses \code{ggplot2} to plot the Moran plot. The plot would be
#' more aesthetically pleasing than the base R version implemented in
#' \code{spdep}. In addition, contours are plotted to show point density on the
#' plot, and the points can be colored by a variable, such as clusters. The
#' contours may also be filled and only influential points plotted. When filled,
#' the viridis E option is used.
#'
#' @inheritParams plotSpatialFeature
#' @inheritParams clusterMoranPlot
#' @param feature Name of one variable to show on the plot. It will be converted
#'   to sentence case on the x axis and lower case in the y axis appended after
#'   "Spatially lagged". One feature at a time since the colors in
#'   \code{color_by} may be specific to this feature (e.g. from
#'   \code{\link{clusterMoranPlot}}).
#' @param color_by Variable to color the points by. It can be the name of a
#'   column in colData, a gene, or the name of a column in the colGeometry
#'   specified in colGeometryName. Or it can be a vector of the same length as
#'   the number of cells/spots in the sample_id of interest.
#' @param plot_singletons Logical, whether to plot items that don't have spatial
#'   neighbors.
#' @param filled Logical, whether to plot filled contours for the
#'   non-influential points and only plot influential points as points.
#' @param colGraphName Name of the \code{colGraph}, the spatial neighborhood
#'   graph used to compute the Moran plot. This is to determine which points are
#'   singletons to plot differently on this plot.
#' @param contour_color Color of the point density contours, which can be
#'   changed so the contours stand out from the points.
#' @param ... Other arguments to pass to \code{\link{geom_density2d}}.
#' @return A ggplot object.
#' @importFrom ggplot2 geom_point aes_string geom_smooth geom_hline geom_vline
#'   geom_density2d scale_shape_manual coord_equal labs geom_density2d_filled
#'   scale_fill_viridis_d scale_x_continuous scale_y_continuous expansion
#'   ggplot_build
#' @importFrom ggplot2 aes
#' @export
moranPlot <- function(sfe, feature, colGraphName, sample_id = NULL,
                      contour_color = "cyan", color_by = NULL,
                      colGeometryName = NULL, annotGeometryName = NULL,
                      plot_singletons = TRUE,
                      filled = FALSE, divergent = FALSE, diverge_center = NULL,
                      name = "MoranPlot", ...) {
  sample_id <- .check_sample_id(sfe, sample_id)
  mp <- .get_feature_metadata(sfe, feature, name, sample_id, colGeometryName,
                              annotGeometryName)[[1]]
  if (isTRUE(is.na(mp))) stop("Moran plot has not been computed for this feature.")
  if (!is.null(color_by)) {
    if (length(color_by) == 1L && is.character(color_by)) {
      # name of something
      if (is.null(annotGeometryName) || !is.null(colGeometryName))
        color_value <- .get_feature_values(sfe, color_by, sample_id,
                                           colGeometryName)
      else {
        ag <- annotGeometry(sfe, annotGeometryName, sample_id)
        color_value <- st_drop_geometry(ag)[ag$sample_id == sample_id,
                                            color_by, drop = FALSE]
      }
    } else if (length(color_by) == sum(colData(sfe)$sample_id == sample_id)) {
      color_value <- color_by
      color_by <- "V1"
    } else {
      stop("color_by must be either the name of a variable in sfe or a vector ",
           "the same length as the number of cells/spots in this sample_id.")
    }
    mp <- cbind(mp, color_value)
  }
  listw <- colGraph(sfe, colGraphName, sample_id)
  is_singleton <- lengths(listw$neighbours) == 0L
  if (filled)
    .moran_ggplot_filled(mp, feature, is_singleton, color_by, plot_singletons,
                         divergent, diverge_center, ...)
  else
    .moran_ggplot(mp, feature, is_singleton, contour_color, color_by, plot_singletons, divergent,
                  diverge_center, ...)
}

#' Plot correlogram
#'
#' Use \code{ggplot2} to plot correlograms computed by
#' \code{\link{runCorrelogram}}, pulling results from \code{rowData}.
#' Correlograms of multiple genes with error bars can be plotted, and they can
#' be colored by any numeric or categorical column in \code{rowData} or a vector
#' with the same length as \code{nrow} of the SFE object. The coloring is useful
#' when the correlograms are clustered to show types of length scales or
#' patterns of decay of spatial autocorrelation. For \code{method = "I"}, the
#' error bars are twice the standard deviation of the estimated Moran's I value.
#'
#' @inheritParams plotSpatialFeature
#' @inheritParams calculateMoransI
#' @inheritParams spdep::sp.correlogram
#' @param color_by Name of a column in \code{rowData(sfe)} or in the
#'   \code{featureData} attribute of \code{colData}, \code{colGeometry}, or
#'   \code{annotGeometry} by which to color the correlogram of each feature.
#' @param plot_signif Logical, whether to plot significance symbols: p < 0.001:
#'   ***, p < 0.01: **, p < 0.05 *, p < 0.1: ., otherwise no symbol. The
#'   p-values are two sided, based on the assumption that the estimated Moran's
#'   I is normally distributed with mean from a randomized version of the data.
#'   The mean and variance come from \code{\link{moran.test}} for Moran's I and
#'   \code{\link{geary.test}} for Geary's C. Take the results with a grain of
#'   salt if the data is not normally distributed.
#' @param p_adj_method Multiple testing correction method as in
#'   \code{\link{p.adjust}}, to correct for multiple testing (number of lags
#'   times number of features) in the Moran's I estimates if \code{plot_signif =
#'   TRUE}.
#' @return A ggplot object.
#' @importFrom ggplot2 theme geom_errorbar element_blank geom_text
#' @importFrom stats p.adjust pnorm symnum
#' @export
plotCorrelogram <- function(sfe, features, sample_id = NULL, method = "I",
                            color_by = NULL,
                            colGeometryName = NULL, annotGeometryName = NULL,
                            plot_signif = TRUE, p_adj_method = "BH",
                            divergent = FALSE, diverge_center = NULL,
                            name = paste("Correlogram", method, sep = "_")) {
  sample_id <- .check_sample_id(sfe, sample_id)
  ress <- .get_feature_metadata(sfe, features, name, sample_id, colGeometryName,
                                annotGeometryName)
  if (!is.null(color_by)) {
    # Different from moranPlot
    if (is.character(color_by) && length(color_by) == 1L) {
      color_value <- .get_feature_metadata(sfe, features, color_by, sample_id,
                                           colGeometryName,
                                           annotGeometryName)
      color_value <- color_value[names(ress)]
    } else if (length(color_by) == length(features)) {
      if (is.null(names(color_by))) names(color_by) <- features
      color_value <- color_by[names(ress)]
    } else {
      stop("color_by must be either the name of a feature in sfe or a vector ",
           "the same length as the number of the features argument.")
    }
  }
  if (method == "corr") {
    dfs <- lapply(seq_along(ress), function(i) {
      res <- ress[[i]]
      if (isTRUE(is.na(res))) return(NA)
      out <- data.frame(lags = seq_along(res),
                        res = res)
      if (length(ress) > 1L) out$feature <- names(ress)[i]
      if (!is.null(color_by)) out$color_by <- color_value[i]
      out
    })
  } else {
    dfs <- lapply(seq_along(ress), function(i) {
      res <- ress[[i]]
      if (isTRUE(is.na(res))) return(NA)
      out <- as.data.frame(res)
      names(out)[names(out) == method] <- "res"
      out$lags <- seq_len(nrow(out))
      out$sd2 <- 2*sqrt(out$variance)
      out$ymin <- out$res - out$sd2
      out$ymax <- out$res + out$sd2
      if (length(features) > 1L) out$feature <- features[[i]]
      if (!is.null(color_by)) out$color_by <- color_value[i]
      out
    })
  }
  is_na_dfs <- vapply(dfs, function(d) isTRUE(is.na(d)), FUN.VALUE = logical(1))
  if (all(is_na_dfs))
    stop("Correlogram has not been computed for any of the features specified ",
         " with method ", method, " for sample ", sample_id)
  if (any(is_na_dfs))
    warning("Correlogram has not been computed for features ",
            paste(features[is_na_dfs], sep = ", "), " with method ", method,
            " for sample ", sample_id)

  dfs <- dfs[!is_na_dfs]
  df <- Reduce(rbind, dfs)
  if (method %in% c("I", "C") && plot_signif) {
    df$z <- (df$res - df$expectation)/sqrt(df$variance)
    df$p <- 2*pnorm(abs(df$z), lower.tail = FALSE)
    df$p_adj <- p.adjust(df$p, method = p_adj_method)
    df$p_symbol <- format(symnum(df$p_adj, corr = FALSE, na = FALSE,
                                 cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                 symbols = c("***", "**", "*", ".", "")))
  }
  lags <- res <- feature <- expectation <- ymin <- ymax <- p_symbol
  if (length(features) > 1L) {
    if (is.null(color_by)) {
      p <- ggplot(df, aes(lags, res, color = feature))
      pal <- .get_pal(df, feature_aes = list(color = "feature"), option = 1,
                      divergent, diverge_center)
      p <- p + pal
      if (method %in% c("I", "C"))
        p <- p + geom_hline(aes(yintercept = expectation, color = feature),
                            linetype = 2, alpha = 0.7)
    }
    else {
      p <- ggplot(df, aes(lags, res, color = color_by, group = feature))
      if (method %in% c("I", "C"))
        p <- p + geom_hline(aes(yintercept = expectation, color = color_by),
                           linetype = 2, alpha = 0.7)
    }
  } else {
    p <- ggplot(df, aes(lags, res))
    if (method %in% c("I", "C"))
      p <- p + geom_hline(aes(yintercept = expectation), linetype = 2, alpha = 0.7)
  }

  p <- p +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = breaks_extended(n = min(max(df$lags), 10), Q = 1:5)) +
    theme(panel.grid.minor.x = element_blank())
  if (method %in% c("I", "C")) {
    p <- p +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2)
    if (plot_signif)
      p <- p + geom_text(aes(y = ymax, label = p_symbol), vjust = 0,
                         show.legend = FALSE)
  } else {
    p <- p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.7)
  }
  if (!is.null(color_by) && length(features) > 1L) {
    pal <- .get_pal(df, feature_aes = list(color = "color_by"), option = 1,
                    divergent, diverge_center)
    p <- p + pal
  }
  p <- p +
    labs(x = "Lags", y = switch(method,
                                corr = "Pearson correlation",
                                I = "Moran's I",
                                C = "Geary's C"))
  p
}

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
ElbowPlot <- function(sce, ndims = 20, reduction = "PCA") {
  # For scater::runPCA results
  percent_var <- attr(reducedDim(sce, reduction), "percentVar")
  if (length(percent_var) < ndims) ndims <- length(percent_var)
  inds <- seq_len(ndims)
  df <- data.frame(PC = inds,
                   pct_var = percent_var[inds])
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
    n2 <- floor(nfeatures/2)
    ord_plus <- order(df$value, decreasing = TRUE)
    ord_minus <- order(df$value, decreasing = FALSE)
    out <- rbind(df[ord_plus[seq_len(n2)],], df[ord_minus[seq_len(n2)],])
  } else {
    ord <- order(abs(df$value), decreasing = TRUE)
    out <- df[ord[seq_len(nfeatures)],]
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
#' @param show_symbol Logical; if the row names of the matrix are Ensembl accessions,
#' indicate whether to show more human readable gene symbols in the plot instead.
#' Ignored if the column specified in \code{symbol_col} is absent from rowData.
#' @param symbol_col If the row names of the gene expression matrix are Ensembl accessions
#' to avoid ambiguity in analysis. If not found in \code{rowData}, then rownames
#' of the gene count matrix will be used.
#' @param balanced Return an equal number of genes with + and - scores. If
#'   FALSE, returns the top genes ranked by the scores absolute values.
#' @param ncol Number of columns in the facetted plot.
#' @return A ggplot object. Loadings for different PCs are plotted in different
#' facets so one ggplot object is returned.
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 facet_wrap scale_y_discrete
#' @importFrom stats reorder
plotDimLoadings <- function(sce, dims = 1:4, nfeatures = 10,
                            show_symbol = TRUE, symbol_col = "symbol",
                            reduction = "PCA",
                            balanced = TRUE, ncol = 2) {
  # For scater::runPCA results
  loadings <- attr(reducedDim(sce, reduction), "rotation")
  is_ensembl <- all(grepl("^ENS", rownames(sce)))
  if (!is_ensembl) show_symbol <- FALSE # I mean, irrelevant
  loading_cols <- paste0("PC", dims)
  df <- cbind(as.data.frame(rowData(sce)[rownames(loadings),]), loadings[,loading_cols])
  if (!symbol_col %in% names(df) || !show_symbol) {
    df$gene_show <- rownames(loadings)
  } else {
    df$gene_show <- df[[symbol_col]]
  }
  df_plt <- lapply(loading_cols, function(p) {
    df_use <- df[,c("gene_show", p)]
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
    facet_wrap(~ PC, scales = "free_y", ncol = 2) +
    scale_y_discrete(labels = function(x) gsub(reg, "", x)) +
    labs(x = "Loading", y = "Gene")
}
