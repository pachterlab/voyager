# 10. What to do with the image when using geom_sf
# Shall I also allow users to plot dimension reductions as features?
# For example, plotting PC1 in space, as opposed to MULTISPATI PC1. I think I'll
# do that, not only for plotting functions, but also for the metrics.
# To do:
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
#' @examples
#' v <- rnorm(10)
#' getDivergeRange(v, diverge_center = 0)
getDivergeRange <- function(values, diverge_center = 0) {
    rg <- range(values, na.rm = TRUE)
    rg_centered <- abs(rg - diverge_center)
    if (!(diverge_center >= rg[1] && diverge_center <= rg[2])) {
        if (diverge_center < rg[1]) {
            pal_begin <- 0.5 + (rg_centered[1]/rg_centered[2])/2
            pal_end <- 1
        } else {
            pal_begin <- 0
            pal_end <- 0.5 - (rg_centered[2]/diverge_center)/2
        }
    } else {
        if (rg_centered[1] < rg_centered[2]) {
            pal_begin <- (rg_centered[2] - rg_centered[1]) / rg_centered[2] / 2
            pal_end <- 1
        } else {
            pal_begin <- 0
            pal_end <- 1 - (rg_centered[1] - rg_centered[2]) / rg_centered[1] / 2
        }
    }
    c(pal_begin, pal_end)
}

.get_applicable <- function(type, fixed) {
    fixed <- .fill_defaults(fixed)
    if (type %in% c("POINT", "MULTIPOINT")) {
        names_use <- c("size", "shape", "alpha", "color")
        if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 0.5
        if (is.na(fixed$color)) fixed$color <- "black"
        shape <- fixed[["shape"]]
        if (!is.null(shape) && shape > 20L) {
            names_use <- c(names_use, "fill")
        }
    } else if (type %in% c("LINESTRING", "MULTILINESTRING")) {
        if (isTRUE(all.equal(0, fixed$size))) fixed$size <- 0.5
        if (is.na(fixed$color)) fixed$color <- "black"
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
        .pal <- switch(option,
            ditto_colors,
            rev(ditto_colors)
        )
        pal_fun <- switch(.aes,
            fill = scale_fill_manual,
            color = scale_color_manual
        )
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
            .pal <- switch(option,
                "roma",
                "bam"
            )
            pal_fun <- switch(.aes,
                fill = scale_fill_scico,
                color = scale_color_scico
            )
            pal <- pal_fun(
                palette = .pal, begin = pal_begin,
                end = pal_end, na.value = "gray"
            )
        } else {
            pal_fun <- switch(.aes,
                fill = scale_fill_distiller,
                color = scale_color_distiller
            )
            .pal <- switch(option,
                "Blues",
                "PuRd"
            )
            pal <- pal_fun(na.value = "gray", palette = .pal, direction = 1)
        }
    }
    pal
}

.fill_defaults <- function(fixed) {
    defaults <- list(
        size = 0, shape = 16,
        linetype = 1, alpha = 1, color = NA, fill = "gray80",
        divergent = FALSE, diverge_center = NULL
    )
    fill <- defaults[setdiff(names(defaults), names(fixed))]
    out <- .drop_null_list(c(fixed, fill))
    if (is.na(out$fill) && "size" %in% names(fill)) {
        out$size <- 0.5
        out$color <- "black"
    }
    out
}

#' @importFrom sf st_drop_geometry st_geometry_type
#' @importFrom ggplot2 ggplot aes_string geom_sf scale_fill_manual
#' scale_color_manual scale_fill_distiller scale_color_distiller geom_polygon
#' geom_segment stat_density2d
#' @importFrom scico scale_fill_scico scale_color_scico
#' @importFrom ggnewscale new_scale_color new_scale_fill
.plot_var_sf <- function(df, annot_df, type, type_annot, feature_aes,
                         feature_fixed, annot_aes, annot_fixed, divergent,
                         diverge_center,annot_divergent, annot_diverge_center,
                         ncol_sample) {
    # Add annotGeometry if present
    if (!is.null(annot_df)) {
        annot_fixed <- .get_applicable(type_annot, annot_fixed)
        annot_fixed <- annot_fixed[setdiff(names(annot_fixed), names(annot_aes))]
        if ("color" %in% names(annot_aes) && annot_fixed$size == 0) {
            annot_fixed$size <- 0.5
        }
        aes_annot <- do.call(aes_string, annot_aes)
        geom_annot <- do.call(geom_sf, c(
            list(mapping = aes_annot, data = annot_df),
            annot_fixed
        ))
        pal_annot <- .get_pal(annot_df, annot_aes, 2, annot_divergent,
                              annot_diverge_center)
    }

    p <- ggplot()
    # Filled polygon annotations go beneath feature plot
    is_annot_filled <- !is.null(annot_df) &&
        ("fill" %in% names(c(annot_aes, annot_fixed))) &&
        type_annot %in% c("POLYGON", "MULTIPOLYGON")
    if ("fill" %in% names(annot_fixed)) {
        is_annot_filled <- is_annot_filled && !is.na(annot_fixed[["fill"]])
    }
    if (is_annot_filled) {
        p <- p + geom_annot
        if (!is.null(pal_annot)) p <- p + pal_annot
    }

    feature_fixed <- .get_applicable(type, feature_fixed)
    feature_fixed <- feature_fixed[setdiff(names(feature_fixed),
                                           names(feature_aes))]

    if ("fill" %in% names(feature_aes) && is_annot_filled) {
        p <- p + new_scale_fill()
    }
    aes_use <- do.call(aes_string, feature_aes)
    geom_use <- do.call(geom_sf, c(
        list(mapping = aes_use, data = df),
        feature_fixed
    ))
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
    if ("sample_id" %in% names(df) && length(unique(df$sample_id)) > 1L) {
        p <- p +
            facet_wrap(~sample_id, ncol = ncol_sample)
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
    # if (type %in% c("POLYGON", "MULTIPOLYGON") || shape >= 20) {
    #  aes_spec <- "fill"
    # }
    if (!.is_discrete(m) && aes_spec %in% c("shape", "linetype")) {
        stop("Shape and linetype are only applicable to discrete variables.")
    }
    aes_spec
}

.get_generalized_geometry_type <- function(g) {
    out <- st_geometry_type(g, by_geometry = FALSE)
    if (out == "GEOMETRY") {
        type_rank <- c(
            "MULTIPOLYGON", "POLYGON", "MULTILINESTRING", "LINESTRING",
            "MULTIPOINT", "POINT", "GEOMETRYCOLLECTION"
        )
        types <- unique(st_geometry_type(g))
        # From most to least general
        out <- type_rank[which(type_rank %in% types)[1]]
    }
    out
}

.wrap_spatial_plots <- function(df, annot_df, type_annot, values, aes_use,
                                annot_aes, annot_fixed, size, shape, linetype,
                                alpha, color, fill, ncol, ncol_sample,
                                divergent, diverge_center, annot_divergent,
                                annot_diverge_center, ...) {
    feature_fixed <- list(
        size = size, shape = shape, linetype = linetype,
        alpha = alpha, color = color, fill = fill
    )
    type <- .get_generalized_geometry_type(df)
    plots <- lapply(names(values), function(n) {
        df[[n]] <- values[[n]]
        feature_aes_name <- .get_feature_aes(df[[n]], type, aes_use, shape)
        feature_aes <- setNames(list(n), feature_aes_name)
        .plot_var_sf(
            df, annot_df, type, type_annot, feature_aes, feature_fixed,
            annot_aes, annot_fixed, divergent, diverge_center,
            annot_divergent, annot_diverge_center, ncol_sample
        )
    })
    if (length(plots) > 1L) {
        out <- wrap_plots(plots, ncol = ncol, ...)
    } else {
        out <- plots[[1]]
    }
    out
}

.plotSpatialFeature <- function(sfe, values, colGeometryName, sample_id, ncol,
                                ncol_sample, annotGeometryName, annot_aes,
                                annot_fixed, aes_use, divergent,
                                diverge_center, annot_divergent,
                                annot_diverge_center, size, shape, linetype,
                                alpha, color, fill, ...) {
    df <- colGeometry(sfe, colGeometryName, sample_id = sample_id)
    if (length(sample_id) > 1L) {
        df$sample_id <- colData(sfe)$sample_id[colData(sfe)$sample_id %in% sample_id]
    }
    # Will use separate ggplots for each feature so each can have its own color scale
    if (!is.null(annotGeometryName)) {
        annot_df <- annotGeometry(sfe, annotGeometryName, sample_id)
        type_annot <- .get_generalized_geometry_type(annot_df)
    } else {
        annot_df <- NULL
        type_annot <- NULL
    }
    .wrap_spatial_plots(
        df, annot_df, type_annot, values, aes_use,
        annot_aes, annot_fixed, size, shape, linetype, alpha,
        color, fill, ncol, ncol_sample, divergent,
        diverge_center, annot_divergent, annot_diverge_center,
        ...
    )
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
#' @inheritParams calculateUnivariate
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
#' @param ncol_sample If plotting multiple samples as facets, how many columns
#'   of such facets. This is distinct from \code{ncols}, which is for multiple
#'   features. When plotting multiple features for multiple samples, then the
#'   result is a multi-panel plot each panel of which is a plot for each feature
#'   facetted by samples.
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
#' @param show_symbol Logical, whether to show human readable gene symbol on the
#'   plot instead of Ensembl IDs when the row names are Ensembl IDs. There must
#'   be a column in \code{rowData(sfe)} called "symbol" for this to work.
#' @param ... Other arguments passed to \code{\link{wrap_plots}}.
#' @return A \code{ggplot2} object if plotting one feature. A \code{patchwork}
#' object if plotting multiple features.
#' @importFrom patchwork wrap_plots
#' @importFrom stats setNames
#' @importMethodsFrom Matrix t
#' @export
#' @examples
#' library(SFEData)
#' library(sf)
#' sfe <- McKellarMuscleData("small")
#' # features can be genes or colData or colGeometry columns
#' plotSpatialFeature(sfe, c("nCounts", rownames(sfe)[1]),
#'     exprs_values = "counts",
#'     colGeometryName = "spotPoly",
#'     annotGeometryName = "tissueBoundary"
#' )
#' # Change fixed aesthetics
#' plotSpatialFeature(sfe, "nCounts",
#'     colGeometryName = "spotPoly",
#'     annotGeometryName = "tissueBoundary",
#'     annot_fixed = list(color = "blue", size = 0.3, fill = NA),
#'     alpha = 0.7
#' )
#' # Make the myofiber segmentations a valid POLYGON geometry
#' ag <- annotGeometry(sfe, "myofiber_simplified")
#' ag <- st_buffer(ag, 0)
#' ag <- ag[!st_is_empty(ag), ]
#' annotGeometry(sfe, "myofiber_simplified") <- ag
#' # Also plot an annotGeometry variable
#' plotSpatialFeature(sfe, "nCounts",
#'     colGeometryName = "spotPoly",
#'     annotGeometryName = "myofiber_simplified",
#'     annot_aes = list(fill = "area")
#' )
plotSpatialFeature <- function(sfe, features, colGeometryName = 1L,
                               sample_id = NULL, ncol = NULL,
                               ncol_sample = NULL,
                               annotGeometryName = NULL,
                               annot_aes = list(), annot_fixed = list(),
                               exprs_values = "logcounts",
                               aes_use = c("fill", "color", "shape", "linetype"),
                               divergent = FALSE, diverge_center = NULL,
                               annot_divergent = FALSE,
                               annot_diverge_center = NULL,
                               size = 0, shape = 16, linetype = 1, alpha = 1,
                               color = NA, fill = "gray80", show_symbol = TRUE,
                               ...) {
    aes_use <- match.arg(aes_use)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    values <- .get_feature_values(sfe, features, sample_id,
        colGeometryName = colGeometryName,
        exprs_values = exprs_values,
        show_symbol = show_symbol
    )
    .plotSpatialFeature(
        sfe, values, colGeometryName, sample_id, ncol,
        ncol_sample, annotGeometryName, annot_aes,
        annot_fixed, aes_use, divergent,
        diverge_center, annot_divergent,
        annot_diverge_center, size, shape, linetype,
        alpha, color, fill, ...
    )
}

.get_graph_df <- function(sfe, MARGIN, sample_id, graph_name, geometry) {
    if (MARGIN == 1L) {
        stop("Not implemented for rowGeometry yet.")
    }
    if (MARGIN == 3L && is.null(geometry)) {
        stop("annotGeometry must be specified.")
    }
    listws <- spatialGraphs(sfe, MARGIN, sample_id = sample_id,
                            name = graph_name)
    sample_inds <- colData(sfe)$sample_id %in% sample_id
    if (is.null(geometry)) {
        coords <- as.data.frame(spatialCoords(sfe)[sample_inds, ])
    } else {
        coords <- as.data.frame(st_coordinates(st_centroid(st_geometry(geometry))))
    }
    if (MARGIN == 2L) {
        coords$sample_id <- colData(sfe)$sample_id[sample_inds]
    } else {
        # sample_id column is required for annotGeometries
        coords$sample_id <- geometry$sample_id
    }
    dfs <- lapply(sample_id, function(s) {
        listw <- listws[[s]]
        cardnb <- card(listw$neighbours)
        n <- length(listw$neighbours)
        neighbors_use <- listw$neighbours[cardnb > 0]
        i <- rep(seq_len(n), cardnb)
        j <- unlist(neighbors_use)
        cu <- coords[coords$sample_id == s, ]
        x_coord <- as.vector(t(cbind(cu[, 1][i], cu[, 1][j])))
        y_coord <- as.vector(t(cbind(cu[, 2][i], cu[, 2][j])))
        df <- data.frame(x = x_coord, y = y_coord)
        df$ID <- rep(seq_len(length(i)), each = 2)
        df <- df2sf(df, geometryType = "LINESTRING")
        df$sample_id <- s
        wts <- listw$weights[cardnb > 0]
        df$weights <- unlist(wts)
        df
    })
    do.call(rbind, dfs)
}

#' @importFrom ggplot2 geom_line
#' @importFrom SpatialFeatureExperiment rowGeometry df2sf
.plot_graph <- function(sfe, MARGIN, sample_id, graph_name, geometry_name,
                        segment_size = 0.5, geometry_size = 0.5, ncol = NULL,
                        weights = FALSE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (!is.null(geometry_name)) {
        gf <- switch(MARGIN,
            rowGeometry,
            colGeometry,
            annotGeometry
        )
        geometry <- gf(sfe, type = geometry_name, sample_id = sample_id)
    } else {
        geometry <- NULL
    }
    df <- .get_graph_df(sfe, MARGIN, sample_id, graph_name, geometry)
    p <- ggplot()
    if (!is.null(geometry_name)) {
        p <- p +
            geom_sf(
                data = geometry, size = geometry_size, fill = NA,
                color = "gray80"
            )
    }
    if (weights) {
        p <- p +
            geom_sf(
                data = df, aes(alpha = weights), size = segment_size,
                color = "black"
            )
    } else {
        p <- p +
            geom_sf(data = df, size = segment_size, color = "black")
    }
    if (length(sample_id) > 1L) {
        p <- p + facet_wrap(~sample_id, ncol = ncol)
    }
    p
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
#' @param weights Whether to plot weights. If \code{TRUE}, then transparency
#'   (alpha) of the segments will represent edge weights.
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SpatialFeatureExperiment spatialGraphs spatialGraph
#' @importFrom spdep card
#' @importFrom sf st_coordinates st_centroid st_geometry
#' @return A ggplot2 object.
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SFEData)
#' library(sf)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' plotColGraph(sfe, colGraphName = "visium", colGeometryName = "spotPoly")
#' # Make the myofiber segmentations a valid POLYGON geometry
#' ag <- annotGeometry(sfe, "myofiber_simplified")
#' ag <- st_buffer(ag, 0)
#' ag <- ag[!st_is_empty(ag), ]
#' annotGeometry(sfe, "myofiber_simplified") <- ag
#' annotGraph(sfe, "myofibers") <-
#'     findSpatialNeighbors(sfe,
#'         type = "myofiber_simplified", MARGIN = 3,
#'         method = "tri2nb", dist_type = "idw"
#'     )
#' plotAnnotGraph(sfe,
#'     annotGraphName = "myofibers",
#'     annotGeometryName = "myofiber_simplified",
#'     weights = TRUE
#' )
plotColGraph <- function(sfe, colGraphName = 1L, colGeometryName = NULL,
                         sample_id = NULL, weights = FALSE, segment_size = 0.5,
                         geometry_size = 0.5, ncol = NULL) {
    .plot_graph(sfe,
        MARGIN = 2L, sample_id = sample_id, graph_name = colGraphName,
        geometry_name = colGeometryName,
        segment_size = segment_size, geometry_size = geometry_size,
        ncol = ncol, weights = weights
    )
}

#' @rdname plotColGraph
#' @export
plotAnnotGraph <- function(sfe, annotGraphName = 1L, annotGeometryName = 1L,
                           sample_id = NULL, weights = FALSE,
                           segment_size = 0.5, geometry_size = 0.5,
                           ncol = NULL) {
    .plot_graph(sfe,
        MARGIN = 3L, sample_id = sample_id, graph_name = annotGraphName,
        geometry_name = annotGeometryName,
        segment_size = segment_size, geometry_size = geometry_size,
        ncol = ncol, weights = weights
    )
}
