.moran_ggplot <- function(mp, feature, is_singleton, contour_color = "cyan",
                          color_by = NULL, plot_singletons = TRUE,
                          divergent = FALSE, diverge_center = NULL, bins = NULL, ...) {
    if (!plot_singletons) {
        mp <- mp[!is_singleton, ]
    }
    x <- wx <- is_inf <- NULL
    if (all(!is_singleton) && plot_singletons) plot_singletons <- FALSE
    p <- ggplot(mp, aes(x = x, y = wx))
    if (plot_singletons) {
        # Need to use listw to check for singletons.
        p <- p +
            geom_point(
                data = mp[is_singleton, ], shape = 21, size = 5, fill = "gray",
                color = "black"
            )
    }
    if (!is.null(color_by)) {
        pal <- .get_pal(mp, list(color = color_by),
            option = 1,
            divergent = divergent, diverge_center = diverge_center
        )
        p <- p + pal
        pts <- geom_point(aes(shape = .data[["is_inf"]],
                              color = .data[[color_by]]))
    } else {
        pts <- geom_point(aes(shape = is_inf), alpha = 0.5)
    }
    p <- p + pts
    # stat_density2d doesn't work when there're too few points
    # Unlikely in real data, but just in case
    # The error doesn't show up until the plot is built.
    p_test <- tryCatch(ggplot_build(p + stat_density2d(bins = bins, ...)),
        error = function(e) {
            warning("Too few points for stat_density2d, not plotting contours.")
        },
        warning = function(w) {
            warning("Too few points for stat_density2d, not plotting contours.")
        }
    )
    if (is(p_test, "ggplot_built")) {
        p <- p + geom_density2d(color = contour_color, bins = bins, ...)
    }
    p <- p +
        geom_smooth(formula = y ~ x, method = "lm") +
        geom_hline(yintercept = mean(mp$wx), lty = 2, color = "gray") +
        geom_vline(xintercept = mean(mp$x), lty = 2, color = "gray") +
        scale_shape_manual(values = c(1, 9)) +
        coord_equal() +
        labs(
            x = feature,
            y = paste("Spatially lagged", feature),
            shape = "Influential"
        )
    p
}

.moran_ggplot_bin2d <- function(mp, feature, is_singleton,
                                plot_singletons = TRUE, bins = 100,
                                binwidth = NULL, hex = FALSE,
                                plot_influential = TRUE) {
    if (!plot_singletons) {
        mp <- mp[!is_singleton, ]
    }
    x <- wx <- is_inf <- NULL
    if (all(!is_singleton) && plot_singletons) plot_singletons <- FALSE
    bin_fun <- if (hex) geom_hex else geom_bin2d
    p <- ggplot(mapping = aes(x = x, y = wx))

    if (plot_influential) {
        p <- p +
            bin_fun(bins = bins, binwidth = binwidth, data = mp[!mp$is_inf, ]) +
            scale_fill_distiller(palette = "Blues", direction = 1) +
            new_scale_fill() +
            bin_fun(bins = bins, binwidth = binwidth, data = mp[mp$is_inf, ]) +
            scale_fill_distiller(palette = "Reds", direction = 1,
                                 name = "count (influential)")
    } else {
        p <- p + bin_fun(bins = bins, binwidth = binwidth, data = mp) +
            scale_fill_distiller(palette = "Blues", direction = 1)
    }

    if (plot_singletons) {
        p <- p +
            geom_point(
                data = mp[is_singleton, ], shape = 21, size = 5, fill = "gray",
                color = "black", alpha = 0.3
            )
    }
    p <- p +
        geom_smooth(formula = y ~ x, method = "lm", data = mp) +
        geom_hline(yintercept = mean(mp$wx), lty = 2, color = "gray") +
        geom_vline(xintercept = mean(mp$x), lty = 2, color = "gray") +
        scale_shape_manual(values = c(1, 9)) +
        coord_equal() +
        labs(
            x = feature,
            y = paste("Spatially lagged", feature)
        )
    p
}

.moran_ggplot_filled <- function(mp, feature, is_singleton, color_by = NULL,
                                 plot_singletons = TRUE, divergent = FALSE,
                                 diverge_center = NULL, bins = NULL, ...) {
    if (!plot_singletons) {
        mp <- mp[!is_singleton, ]
    }
    x <- wx <- is_inf <- NULL
    p <- ggplot(mp, aes(x = x, y = wx)) +
        geom_density2d_filled(show.legend = FALSE, bins = bins, ...)
    if (plot_singletons) {
        p <- p +
            geom_point(
                data = mp[is_singleton & mp$is_inf, ], shape = 21, size = 5,
                fill = "blue", color = "cornflowerblue"
            )
    }
    mp_inf <- mp[mp$is_inf, ]
    if (!is.null(color_by)) {
        pal <- .get_pal(mp_inf, list(color = color_by),
            option = 1,
            divergent = divergent, diverge_center = diverge_center
        )
        p <- p + pal
        pts <- geom_point(
            data = mp_inf, aes(color = .data[[color_by]]),
            shape = 9
        )
    } else {
        pts <- geom_point(data = mp_inf, shape = 9, color = "cornflowerblue")
    }
    p +
        geom_smooth(formula = y ~ x, method = "lm") +
        geom_hline(yintercept = mean(mp$wx), lty = 2, color = "cornflowerblue") +
        geom_vline(xintercept = mean(mp$x), lty = 2, color = "cornflowerblue") +
        pts +
        scale_fill_viridis_d(option = "E") +
        coord_equal() +
        scale_x_continuous(expand = expansion()) +
        scale_y_continuous(expand = expansion()) +
        labs(
            x = feature,
            y = paste("Spatially lagged", feature),
            shape = "Influential"
        ) +
        theme(panel.background = element_rect(fill = scales::viridis_pal(option = "E")(255)[1]),
              panel.grid = element_blank())
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
#' @inheritParams plotCorrelogram
#' @inheritParams plotCellBin2D
#' @inheritParams plotDimLoadings
#' @param sample_id One sample_id for the sample whose graph to plot.
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
#' @param binned Logical, whether to plot 2D histograms. This argument has
#'   precedence to \code{filled}.
#' @param graphName Name of the \code{colGraph} or \code{annotGraph}, the
#'   spatial neighborhood graph used to compute the Moran plot. This is to
#'   determine which points are singletons to plot differently on this plot.
#' @param contour_color Color of the point density contours, which can be
#'   changed so the contours stand out from the points.
#' @param plot_influential Logical, whether to plot influential points with
#'   different palette if \code{binned = TRUE}.
#' @param name Name under which the Moran plot results are stored. By default
#'   "moran.plot".
#' @param bins_contour Number of bins in the point density contour. Use a
#'   smaller number to make sparser contours.
#' @param ... Other arguments to pass to \code{\link{geom_density2d}}.
#' @return A ggplot object.
#' @importFrom ggplot2 geom_point geom_smooth geom_hline geom_vline
#'   geom_density2d scale_shape_manual coord_equal labs geom_density2d_filled
#'   scale_fill_viridis_d scale_x_continuous scale_y_continuous expansion
#'   ggplot_build aes geom_hex geom_bin2d
#' @importFrom SpatialFeatureExperiment localResult
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' library(bluster)
#' library(scater)
#' sfe <- McKellarMuscleData("full")
#' sfe <- sfe[, colData(sfe)$in_tissue]
#' sfe <- logNormCounts(sfe)
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' sfe <- runUnivariate(sfe, type = "moran.plot", features = "Myh1",
#'                      swap_rownames = "symbol")
#' clust <- clusterMoranPlot(sfe, "Myh1", BLUSPARAM = KmeansParam(2),
#'                           swap_rownames = "symbol")
#' moranPlot(sfe, "Myh1", graphName = "visium", color_by = clust[, 1],
#'           swap_rownames = "symbol")
moranPlot <- function(sfe, feature, graphName = 1L, sample_id = "all",
                      contour_color = "cyan", color_by = NULL,
                      colGeometryName = NULL, annotGeometryName = NULL,
                      plot_singletons = TRUE, binned = FALSE,
                      filled = FALSE, divergent = FALSE, diverge_center = NULL,
                      swap_rownames = NULL,
                      bins = 100, binwidth = NULL, hex = FALSE,
                      plot_influential = TRUE, bins_contour = NULL,
                      name = "moran.plot", ...) {
    l <- .deprecate_show_symbol("moranPlot", show_symbol, swap_rownames)
    show_symbol <- l[[1]]; swap_rownames <- l[[2]]

    sample_id <- .check_sample_id(sfe, sample_id)
    # Change as moran.plot has been moved to localResults.
    not_geometry <- is.null(colGeometryName) && is.null(annotGeometryName)
    if (not_geometry) feature <- .symbol2id(sfe, feature, swap_rownames)
    mp <- localResult(sfe,
        type = name, feature = feature,
        sample_id = sample_id, colGeometryName = colGeometryName,
        annotGeometryName = annotGeometryName
    )
    if (!is.null(swap_rownames) && not_geometry) {
        if (feature %in% rownames(sfe) && swap_rownames %in% colnames(rowData(sfe))) {
            feature <- rowData(sfe)[feature, swap_rownames]
        }
    }

    if (isTRUE(is.na(mp))) stop("Moran plot has not been computed for this feature.")
    use_col <- is.null(annotGeometryName) || !is.null(colGeometryName)
    mar <- if (use_col) 2L else 3L
    if (use_col) {
        length_ref <- sum(colData(sfe)$sample_id == sample_id)
    } else {
        ag <- annotGeometry(sfe, annotGeometryName, sample_id)
        ag <- .rm_empty_geometries(ag, 3)
        length_ref <- nrow(ag)
    }
    if (!is.null(color_by)) {
        if (length(color_by) == 1L && is.character(color_by)) {
            # name of something
            if (use_col) {
                color_value <- .get_feature_values(
                    sfe, color_by, sample_id,
                    colGeometryName
                )
            } else {
                color_value <- st_drop_geometry(ag)[ag$sample_id == sample_id,
                    color_by,
                    drop = FALSE
                ]
            }
        } else if (length(color_by) == length_ref) {
            color_value <- color_by
            color_by <- "color_value"
        } else {
            stop(
                "color_by must be either the name of a variable in sfe or a vector ",
                "the same length as the number of cells/spots in this sample_id."
            )
        }
        mp <- cbind(mp, color_value)
    }
    listw <- spatialGraph(sfe, type = graphName, MARGIN = mar,
                          sample_id = sample_id)
    is_singleton <- vapply(listw$neighbours, min, FUN.VALUE = integer(1)) == 0L
    if (binned) {
        .moran_ggplot_bin2d(mp, feature, is_singleton, plot_singletons,
                            bins = bins, binwidth = binwidth, hex = hex,
                            plot_influential = plot_influential)
    } else if (filled) {
        .moran_ggplot_filled(
            mp, feature, is_singleton, color_by, plot_singletons,
            divergent, diverge_center, bins = bins_contour, ...
        )
    } else {
        .moran_ggplot(
            mp, feature, is_singleton, contour_color, color_by, plot_singletons,
            divergent, diverge_center, bins = bins_contour, ...
        )
    }
}

.get_color_by <- function(sfe, features, color_by, sample_id,
                          colGeometryName, annotGeometryName, reducedDimName,
                          show_symbol, swap_rownames) {
    if (is.data.frame(color_by)) return(NULL)
    if (!is.null(color_by)) {
        # Different from moranPlot
        if (is.character(color_by) && length(color_by) == 1L) {
            color_value <- .get_feature_metadata(
                sfe, features, color_by, sample_id,
                colGeometryName, annotGeometryName, reducedDimName,
                show_symbol, swap_rownames
            )
            color_value <- color_value[features]
        } else if (length(color_by) == length(features)) {
            if (is.null(names(color_by))) names(color_by) <- features
            color_value <- color_by[features]
        } else {
            stop(
                "color_by must be either the name of a feature in sfe or a vector ",
                "the same length as the features argument."
            )
        }
        color_value
    } else NULL
}

.get_plot_correlogram_df <- function(sfe, features, sample_id, method, color_by,
                                     colGeometryName, annotGeometryName,
                                     reducedDimName, name,
                                     show_symbol, swap_rownames) {
    ress <- .get_feature_metadata(
        sfe, features, name, sample_id, colGeometryName,
        annotGeometryName, reducedDimName, show_symbol, swap_rownames
    )
    color_value <- .get_color_by(sfe, features, color_by, sample_id,
                              colGeometryName, annotGeometryName, reducedDimName,
                              show_symbol, swap_rownames)
    if (method == "corr") {
        dfs <- lapply(seq_along(ress), function(i) {
            res <- ress[[i]]
            if (isTRUE(is.na(res))) {
                return(NA)
            }
            out <- data.frame(
                lags = seq_along(res),
                res = res
            )
            if (length(ress) > 1L) out$feature <- names(ress)[i]
            if (!is.null(color_value)) out$color_by <- color_value[i]
            out
        })
    } else {
        dfs <- lapply(seq_along(ress), function(i) {
            res <- ress[[i]]
            if (isTRUE(is.na(res))) {
                return(NA)
            }
            out <- as.data.frame(res)
            names(out)[names(out) == method] <- "res"
            out$lags <- seq_len(nrow(out))
            out$sd2 <- 2 * sqrt(out$variance)
            out$ymin <- out$res - out$sd2
            out$ymax <- out$res + out$sd2
            out$feature <- names(ress)[[i]]
            if (!is.null(color_value)) out$color_by <- color_value[i]
            out
        })
    }
    is_na_dfs <- vapply(dfs, function(d) isTRUE(is.na(d)),
                        FUN.VALUE = logical(1))
    if (all(is_na_dfs)) {
        stop(
            "Correlogram has not been computed for any of the features specified ",
            " with method ", method, " for sample ", sample_id
        )
    }
    if (any(is_na_dfs)) {
        warning(
            "Correlogram has not been computed for features ",
            paste(features[is_na_dfs], sep = ", "), " with method ", method,
            " for sample ", sample_id
        )
    }

    dfs <- dfs[!is_na_dfs]
    do.call(rbind, dfs)
}

.dimred_feature_order <- function(df) {
    # Make sure that PC2 follows PC1, not PC10
    feature_num <- sort(as.numeric(unique(gsub("^[A-Z]*", "", df$feature))))
    pre <- unique(gsub("\\d+$", "", df$feature))
    df$feature <- factor(df$feature, levels = paste0(pre, feature_num))
    df
}
#' Plot correlogram
#'
#' Use \code{ggplot2} to plot correlograms computed by
#' \code{\link{runUnivariate}}, pulling results from \code{rowData}.
#' Correlograms of multiple genes with error bars can be plotted, and they can
#' be colored by any numeric or categorical column in \code{rowData} or a vector
#' with the same length as \code{nrow} of the SFE object. The coloring is useful
#' when the correlograms are clustered to show types of length scales or
#' patterns of decay of spatial autocorrelation. For \code{method = "I"}, the
#' error bars are twice the standard deviation of the estimated Moran's I value.
#'
#' @inheritParams plotSpatialFeature
#' @inheritParams calculateUnivariate
#' @inheritParams spdep::sp.correlogram
#' @inheritParams plotDimLoadings
#' @inheritParams getParams
#' @param features Features to plot, must be in rownames of the gene count
#'   matrix, colnames of colData or a colGeometry, colnames of cell embeddings
#'   in \code{reducedDim}, or numeric indices of dimension reduction components.
#' @param color_by Name of a column in \code{rowData(sfe)} or in the
#'   \code{featureData} of \code{colData} (see \code{\link{colFeatureData}}),
#'   \code{colGeometry}, or \code{annotGeometry} by which to color the
#'   correlogram of each feature. Alternatively, a vector of the same length as
#'   \code{features}, or a data frame from \code{\link{clusterCorrelograms}}.
#' @param colGeometryName Name of a \code{colGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute the metric. Use
#'   \code{\link{colGeometryNames}} to look up names of the \code{sf} data
#'   frames associated with cells/spots.
#' @param plot_signif Logical, whether to plot significance symbols: p < 0.001:
#'   ***, p < 0.01: **, p < 0.05 *, p < 0.1: ., otherwise no symbol. The
#'   p-values are two sided, based on the assumption that the estimated Moran's
#'   I is normally distributed with mean from a randomized version of the data.
#'   The mean and variance come from \code{\link{moran.test}} for Moran's I and
#'   \code{\link{geary.test}} for Geary's C. Take the results with a grain of
#'   salt if the data is not normally distributed.
#' @param facet_by Whether to facet by sample_id (default) or features. If
#'   facetting by sample_id, then different features will be plotted in the same
#'   facet for comparison. If facetting by features, then different samples will
#'   be compared for each feature. Ignored if only one sample is specified.
#' @param ncol Number of columns if facetting.
#' @param p_adj_method Multiple testing correction method as in
#'   \code{\link{p.adjust}}, to correct for multiple testing (number of lags
#'   times number of features) in the Moran's I estimates if \code{plot_signif =
#'   TRUE}.
#' @param name Name under which the correlogram results are stored, which is by
#' default "sp.correlogram".
#' @return A ggplot object.
#' @importFrom ggplot2 theme geom_errorbar element_blank geom_text
#' @importFrom stats p.adjust pnorm symnum
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SFEData)
#' library(bluster)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- logNormCounts(sfe)
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' inds <- c(1, 3, 4, 5)
#' features <- rownames(sfe)[inds]
#' sfe <- runUnivariate(sfe,
#'     type = "sp.correlogram", features = features,
#'     exprs_values = "counts", order = 5
#' )
#' clust <- clusterCorrelograms(sfe,
#'     features = features,
#'     BLUSPARAM = KmeansParam(2)
#' )
#' # Color by features
#' plotCorrelogram(sfe, features)
#' # Color by something else
#' plotCorrelogram(sfe, features, color_by = clust$cluster)
#' # Facet by features
#' plotCorrelogram(sfe, features, facet_by = "features")
plotCorrelogram <- function(sfe, features, sample_id = "all", method = "I",
                            color_by = NULL,
                            facet_by = c("sample_id", "features"),
                            ncol = NULL,
                            colGeometryName = NULL, annotGeometryName = NULL,
                            reducedDimName = NULL,
                            plot_signif = TRUE, p_adj_method = "BH",
                            divergent = FALSE, diverge_center = NULL,
                            swap_rownames = NULL,
                            name = "sp.correlogram") {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (length(sample_id) > 1L || length(features) > 1L) {
        facet_by <- match.arg(facet_by)
    }
    facet_sample <- length(sample_id) > 1L && facet_by == "sample_id"
    facet_feature <- length(features) > 1L && facet_by == "features"
    name <- paste(name, method, sep = "_")
    is_dimred <- is.null(colGeometryName) && is.null(annotGeometryName) &&
        !is.null(reducedDimName) && length(features) > 1L
    df <- lapply(sample_id, function(s) {
        o <- .get_plot_correlogram_df(
            sfe, features, s, method, color_by,
            colGeometryName, annotGeometryName, reducedDimName, name,
            !is.null(swap_rownames), swap_rownames
        )
        o$sample_id <- s
        o
    })
    if (length(sample_id) > 1L) {
        df <- do.call(rbind, df)
    } else {
        df <- df[[1]]
    }
    if (is_dimred){
        df <- .dimred_feature_order(df)
    }
    if (is.data.frame(color_by)) {
        if (length(unique(color_by$cluster)) < 2L) color_by <- NULL
        else {
            df <- merge(df, color_by, by = c("feature", "sample_id"), all.x = TRUE)
            names(df)[names(df) == "cluster"] <- "color_by"
        }
    }
    if (method %in% c("I", "C") && plot_signif) {
        df$z <- (df$res - df$expectation) / sqrt(df$variance)
        df$p <- 2 * pnorm(abs(df$z), lower.tail = FALSE)
        df$p_adj <- p.adjust(df$p, method = p_adj_method)
        df$p_symbol <- format(symnum(df$p_adj,
            corr = FALSE, na = FALSE,
            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", ".", "")
        ))
    }
    lags <- res <- feature <- expectation <- ymin <- ymax <- p_symbol <- NULL
    group_sample <- !facet_sample && length(sample_id) > 1L
    if (length(features) > 1L) {
        if (is.null(color_by)) {
            if (group_sample) {
                p <- ggplot(df, aes(lags, res, color = sample_id))
                pal <- .get_pal(df,
                    feature_aes = list(color = "sample_id"), option = 1,
                    divergent, diverge_center
                )
            } else {
                p <- ggplot(df, aes(lags, res, color = feature))
                if (is_dimred) {
                    pal <- scale_color_viridis_d(end = 0.9, option = "E")
                } else {
                    pal <- .get_pal(df,
                                    feature_aes = list(color = "feature"), option = 1,
                                    divergent, diverge_center
                    )
                }
            }
            p <- p + pal
            if (method %in% c("I", "C")) {
                if (group_sample) {
                    p <- p + geom_hline(aes(yintercept = expectation,
                                            color = sample_id),
                        linetype = 2, alpha = 0.7
                    )
                } else {
                    p <- p + geom_hline(aes(yintercept = expectation,
                                            color = feature),
                        linetype = 2, alpha = 0.7
                    )
                }
            }
        } else {
            if (group_sample) {
                p <- ggplot(df, aes(lags, res, color = color_by,
                                    linetype = sample_id))
            } else {
                p <- ggplot(df, aes(lags, res, color = color_by,
                                    group = feature))
            }
            if (method %in% c("I", "C")) {
                p <- p + geom_hline(aes(yintercept = expectation,
                                        color = color_by),
                    linetype = 2, alpha = 0.7
                )
            }
        }
    } else {
        if (group_sample) {
            p <- ggplot(df, aes(lags, res, color = sample_id)) +
                .get_pal(df,
                    feature_aes = list(color = "sample_id"), option = 1,
                    divergent, diverge_center
                )
        } else {
            p <- ggplot(df, aes(lags, res))
        }
        if (method %in% c("I", "C")) {
            if (group_sample) {
                p <- p + geom_hline(aes(yintercept = expectation,
                                        color = sample_id),
                    linetype = 2, alpha = 0.7
                )
            } else {
                p <- p + geom_hline(aes(yintercept = expectation), linetype = 2,
                                    alpha = 0.7)
            }
        }
    }
    p <- p +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = breaks_extended(n = min(max(df$lags), 10),
                                                    Q = 1:5)) +
        theme(panel.grid.minor.x = element_blank())
    if (method %in% c("I", "C")) {
        p <- p +
            geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2)
        if (plot_signif) {
            p <- p + geom_text(aes(y = ymax, label = p_symbol),
                vjust = 0,
                show.legend = FALSE
            )
        }
    } else {
        p <- p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.7)
    }
    if (!is.null(color_by) && length(features) > 1L) {
        name_use <- if (is.data.frame(color_by)) "cluster" else "color_by"
        pal <- .get_pal(df,
                        feature_aes = list(color = "color_by"), option = 1,
                        divergent, diverge_center, name = name_use
        )
        p <- p + pal
    }
    p <- p +
        labs(x = "Lags", y = switch(method,
            corr = "Pearson correlation",
            I = "Moran's I",
            C = "Geary's C"
        ))
    if (facet_feature) {
        p <- p + facet_wrap(~feature, ncol = ncol)
    }
    if (facet_sample) {
        p <- p + facet_wrap(~sample_id, ncol = ncol)
    }
    p
}

.get_plot_mc_df <- function(sfe, features, sample_id, name,
                            colGeometryName, annotGeometryName, reducedDimName,
                            show_symbol, swap_rownames) {
    # Ah, the weight of tradition. .get_feature_metadata only works for one sample at a time
    # As a result, this function deals with one sample at a time.
    ress <- .get_feature_metadata(sfe, features,
        name = paste0(name, "_res"),
        sample_id = sample_id, colGeometryName,
        annotGeometryName, reducedDimName,
        show_symbol, swap_rownames
    )
    res_stats <- .get_feature_metadata(sfe, features,
        name = paste0(name, "_statistic"),
        sample_id = sample_id, colGeometryName,
        annotGeometryName, reducedDimName,
        show_symbol, swap_rownames
    )
    dfs <- lapply(seq_along(ress), function(i) {
        if (isTRUE(is.na(ress[[i]]))) {
            return(NA)
        }
        res_use <- ress[[i]]
        res_use <- res_use[-length(res_use)]
        data.frame(
            res = res_use,
            statistic = res_stats[[i]],
            feature = names(ress)[i],
            sample_id = sample_id
        )
    })
    dfs <- dfs[!.is_na_list(dfs)]
    if (!length(dfs)) {
        stop("None of the features have the specified MC computed.")
    }
    df <- do.call(rbind, dfs)
    if (is.null(colGeometryName) && is.null(annotGeometryName) &&
        !is.null(reducedDimName) && length(features) > 1L){
        df <- .dimred_feature_order(df)
    }
    df
}

#' Plot Moran/Geary Monte Carlo results
#'
#' Plot the simulations as a density plot or histogram compared to the observed
#' Moran's I or Geary's C, with ggplot2 so it looks nicer. Unlike the plotting
#' function in \code{spdep}, this function can also plot the same feature in
#' different samples as facets or plot different features or samples together
#' for comparison.
#'
#' @inheritParams plotCorrelogram
#' @param ptype Plot type, one of "density", "histogram", or "freqpoly".
#' @param name Name under which the Monte Carlo results are stored, which
#'   defaults to "moran.mc". For Geary's C Monte Carlo, the default is
#'   "geary.mc".
#' @param ... Other arguments passed to \code{\link{geom_density}},
#'   \code{\link{geom_histogram}}, or \code{\link{geom_freqpoly}}, depending on
#'   \code{ptype}.
#' @return A \code{ggplot2} object.
#' @importFrom ggplot2 geom_density geom_histogram geom_freqpoly
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' sfe <- colDataUnivariate(sfe, type = "moran.mc", "nCounts", nsim = 100)
#' plotMoranMC(sfe, "nCounts")
plotMoranMC <- function(sfe, features, sample_id = "all",
                        facet_by = c("sample_id", "features"), ncol = NULL,
                        colGeometryName = NULL, annotGeometryName = NULL,
                        reducedDimName = NULL,
                        ptype = c("density", "histogram", "freqpoly"),
                        swap_rownames = NULL,
                        name = "moran.mc", ...) {
    ptype <- match.arg(ptype)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (length(sample_id) > 1L || length(features) > 1L) {
        facet_by <- match.arg(facet_by)
    }
    facet_sample <- length(sample_id) > 1L && facet_by == "sample_id"
    facet_feature <- length(features) > 1L && facet_by == "features"
    group_sample <- !facet_sample && length(sample_id) > 1L

    dens_geom <- switch(ptype,
        density = geom_density,
        histogram = geom_histogram,
        freqpoly = geom_freqpoly
    )

    df <- lapply(sample_id, function(s) {
        .get_plot_mc_df(
            sfe, features, s, name,
            colGeometryName,
            annotGeometryName, reducedDimName,
            !is.null(swap_rownames), swap_rownames
        )
    })
    if (length(sample_id) > 1L) {
        df <- do.call(rbind, df)
    } else {
        df <- df[[1]]
    }
    p <- ggplot(df)
    statistic <- feature <- res <- NULL
    if ((length(sample_id) == 1L && (length(features) == 1L || facet_feature)) ||
        (length(features) == 1L && facet_sample)) {
        p <- ggplot(df) +
            geom_vline(aes(xintercept = statistic))
    }
    if (length(features) > 1L && (length(sample) == 1L || facet_sample)) {
        if (ptype == "histogram") {
            stop("Histograms are not supported when multiple colors are used.")
        }
        p <- ggplot(df, aes(color = feature)) +
            geom_vline(aes(xintercept = statistic, color = feature))
    }
    if (length(sample_id) > 1L && (length(features) == 1L || facet_feature)) {
        if (ptype == "histogram") {
            stop("Histograms are not supported when multiple colors are used.")
        }
        p <- ggplot(df, aes(color = sample_id)) +
            geom_vline(aes(xintercept = statistic, color = sample_id))
    }
    method_show <- if (grepl("[mM]oran", name)) "Moran's I" else "Geary's C"
    if (is.null(colGeometryName) && is.null(annotGeometryName) &&
        !is.null(reducedDimName)) {
        pal <- scale_color_viridis_d(end = 0.9, option = "E")
    } else {
        pal <- scale_color_manual(values = ditto_colors)
    }
    p <- p + dens_geom(aes(res), ...) + pal +
        labs(x = paste0("Monte-Carlo simulation of ", method_show))
    if (facet_feature) p <- p + facet_wrap(~feature)
    if (facet_sample) p <- p + facet_wrap(~sample_id)
    p
}
