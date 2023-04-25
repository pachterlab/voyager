#' Compute variograms
#'
#' Wrapper of \code{automap::autofitVariogram} to facilitate computing
#' variograms for multiple genes in SFE objects as an EDA tool. These functions
#' are written to conform to a uniform format for univariate methods to be
#' called internally. These functions are not exported, but the documentation is
#' written to show users the extra arguments to use when alling
#' \code{calculateUnivariate} or \code{runUnivariate}.
#'
#' @inheritParams gstat::variogram
#' @param x A numeric vector whose variogram is computed.
#' @param y For bivariate, another numeric vector whose variogram is computed.
#' @param coords_df A \code{sf} data frame with the geometry and regressors for
#'   variogram modeling.
#' @param formula A formula defining the response vector and (possible)
#'   regressors, in case of absence of regressors, use x ~ 1.
#' @param scale Logical, whether to scale \code{x}. Defaults to \code{TRUE} so
#'   the variogram is easier to interpret and is more comparable between
#'   features with different magnitudes when the length scale of spatial
#'   autocorrelation is of interest.
#' @param ... Other arguments passed to \code{automap::autofitVariogram} such as
#'   \code{model} and \code{\link{variogram}} such as \code{alpha} for
#'   anisotropy. Note that \code{gstat} does not fit ansotropic models and you
#'   will get a warning if you specify \code{alpha}. Nevertheless, plotting the
#'   empirical anisotropic variograms and comparing them to the variogram fitted
#'   to the entire dataset can be a useful EDA tool.
#' @return An \code{autofitVariogram} object.
#' @name variogram-internal
.variogram <- function(x, coords_df, formula = x ~ 1, scale = TRUE, ...) {
    coords_df$x <- x
    if (scale) coords_df$x <- scale(coords_df$x)
    dots <- list(...)
    # Deal with alpha myself and fit a global variogram to avoid further gstat warnings
    have_alpha <- "alpha" %in% names(dots)
    if (have_alpha) {
        empirical <- gstat::variogram(formula, data = coords_df, alpha = dots$alpha)
        dots$alpha <- NULL
    }
    out <- do.call(automap::autofitVariogram,
                   c(list(formula = formula, input_data = coords_df,
                          map = FALSE, cloud = FALSE), dots))
    if (have_alpha) {
        out$exp_var <- empirical
    }
    out
}

#' @rdname variogram-internal
.variogram_bv <- function(x, y, coords_df, scale = TRUE, map = FALSE, ...) {
    # should take x as a matrix. All pairwise cross variograms for columns of x are computed.
    # To be consistent with other functions, x has genes as rows
    # should already be reformatted in calculateBivariate
    dots <- list(...)
    use_formula <- "formula" %in% names(dots)
    if (use_formula) {
        form_str <- deparse(dots[["formula"]])
        form_str <- gsub("\\s+",, "", form_str)
        form_split <- strsplit(form_str, "~")[[1]]
        lhs <- form_split[1]
        if (lhs != "") {
            message("Left hand side of formula is ignored.")
        }
        dots$formula <- NULL
    }
    if (!is.null(y)) {
        do_message <- nrow(x) > 1L || nrow(y) > 1L
        x <- rbind(x, y)
        rownames(x) <- make.names(rownames(x)) # in case of duplicates
        if (do_message)
            message("Cross correlograms within columns of x and within columns of y are also computed.")
    }
    ns <- rownames(x)
    x <- t(x)
    if (scale) x <- scale(x)
    coords_df <- st_as_sf(cbind(coords_df, x))

    g <- NULL
    for (n in ns) {
        if (use_formula) {
            formula_use <- paste(c(n, "~", form_split[-1]), collapse = " ") |>
                as.formula()
        } else {
            formula_use <- as.formula(paste0(n, " ~ 1"))
        }
        g <- gstat::gstat(g, id = n, formula = formula_use, data = coords_df)
    }
    do.call(gstat::variogram, c(list(object = g, map = map), dots))
}

.vbv_fun <- function(map) {
    if (!map) {
        function(x, y, coords_df, scale = TRUE, ...)
            .variogram_bv(x, y, coords_df, scale, map = FALSE, ...)
    }
    else {
        function(x, y, coords_df, width, cutoff, scale = TRUE, ...)
            .variogram_bv(x, y, coords_df, scale, map = TRUE, width = width, cutoff = cutoff, ...)
    }
}

#' @rdname variogram-internal
.cross_variogram <- .vbv_fun(map = FALSE)

#' @rdname variogram-internal
.cross_variogram_map <- .vbv_fun(map = TRUE)

#' @rdname variogram-internal
.variogram_map <- function(x, coords_df, formula = x ~ 1, width, cutoff, scale = TRUE, ...) {
    coords_df$x <- x
    if (scale) coords_df$x <- scale(coords_df$x)
    gstat::variogram(formula, data = coords_df, width = width, cutoff = cutoff,
                     map = TRUE, ...)
}

variogram <- SFEMethod(package = "automap", variate = "uni", scope = "global",
                       default_attr = NA, name = "variogram", title = "Variogram",
                       fun = .variogram,
                       reorganize_fun = .other2df,
                       use_graph = FALSE)

cross_variogram <- SFEMethod(package = "gstat", variate = "bi", scope = "global",
                             default_attr = NA, name = "cross_variogram",
                             title = "Cross variogram",
                             fun = .cross_variogram,
                             reorganize_fun = function(out, name, ...) out,
                             use_graph = FALSE, use_matrix = TRUE)

cross_variogram_map <- SFEMethod(
    package = "gstat", variate = "bi", scope = "global",
    default_attr = NA, name = "cross_variogram_map", title = "Cross variogram map",
    fun = .cross_variogram_map,
    reorganize_fun = function(out, name, ...) as.data.frame(out$map),
    use_graph = FALSE, use_matrix = TRUE
)

variogram_map <- SFEMethod(package = "gstat", variate = "uni", scope = "global",
                           default_attr = NA, name = "variogram_map",
                           title = "Variogram map",
                           fun = .variogram_map,
                           reorganize_fun = .variogram_map2df,
                           use_graph = FALSE)

.get_plot_variogram_df <- function(sfe, features, sample_id, color_by,
                                   colGeometryName, annotGeometryName,
                                   reducedDimName, name,
                                   show_symbol, swap_rownames) {
    # For one sample
    ress <- .get_feature_metadata(
        sfe, features, name, sample_id, colGeometryName,
        annotGeometryName, reducedDimName, show_symbol, swap_rownames
    )
    color_value <- .get_color_by(sfe, features, color_by, sample_id,
                                 colGeometryName, annotGeometryName, reducedDimName,
                                 show_symbol, swap_rownames)
    exp_vars <- lapply(seq_along(ress), function(i) {
        df <- ress[[i]]$exp_var
        df$feature <- names(ress)[i]
        if (!is.null(color_value)) df$color_by <- color_value[i]
        df
    })
    exp_vars <- do.call(rbind, exp_vars)
    exp_vars$sample_id <- sample_id

    var_models <- lapply(seq_along(ress), function(i) {
        df <- gstat::variogramLine(ress[[i]]$var_model,
                                   maxdist = max(ress[[i]]$exp_var$dist))
        df$feature <- names(ress)[i]
        if (!is.null(color_value)) df$color_by <- color_value[i]
        df
    })
    var_models <- do.call(rbind, var_models)
    var_models$sample_id <- sample_id

    # For text label on model, nugget, sill, and range
    model_labels <- lapply(seq_along(ress), function(i) {
        df <- ress[[i]]$var_model
        data.frame(
            label = paste0("Model: ", df$model[2],
                           "\nNugget: ", format(df$psill[1], digits = 3),
                           "\nSill: ", format(sum(df$psill), digits = 3),
                           "\nRange: ", format(df$range[2], digits = 3)),
            feature = names(ress)[i]
        )
    })
    model_labels <- do.call(rbind, model_labels)
    model_labels$sample_id <- sample_id
    list(exp_vars = exp_vars,
         var_models = var_models,
         model_labels = model_labels)
}

#' Plot variogram
#'
#' This function plots the variogram of a feature and its fitted variogram
#' models, showing the nugget, range, and sill of the model. Unlike the plotting
#' functions in package \code{automap} that uses \code{lattice}, this function
#' uses \code{ggplot2} to make prettier and more customizable plots.
#'
#' @inheritParams plotCorrelogram
#' @param group Which of samples, features, and angles to show in the same facet
#'   for comparison when there are multiple. Default to "none", meaning each
#'   facet will contain one variogram. When grouping multiple variograms in the
#'   same facet, the text with model, nugget, sill, and range of the variograms
#'   will not be shown.
#' @param use_lty Logical, whether to use linetype or point shape to distinguish
#'   between the different features or samples in the same facet. If
#'   \code{FALSE}, then the different features or samples are not distinguished
#'   and the patterns are shown only.
#' @param show_np Logical, whether to show number of pairs of cells at each
#'   distance bin.
#' @return A \code{ggplot} object. The empirical variogram at each distance bin
#'   is plotted as points, and the fitted variogram model is plotted as a line
#'   for each feature. The number next to each point is the number of pairs of
#'   cells in that distance bin.
#' @seealso plotVariogramMap
#' @export
#' @importFrom rlang !! sym
#' @importFrom utils modifyList
#' @importFrom ggplot2 scale_color_viridis_d
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData()
#' sfe <- colDataUnivariate(sfe, "variogram", features = "nCounts", model = "Sph")
#' plotVariogram(sfe, "nCounts")
#' # Anisotropy, will get a message
#' sfe <- colDataUnivariate(sfe, "variogram", features = "nCounts",
#' model = "Sph", alpha = c(30, 90, 150), name = "variogram_anis")
#' # Facet by angles by default
#' plotVariogram(sfe, "nCounts", name = "variogram_anis")
#' # Plot angles with different colors
#' plotVariogram(sfe, "nCounts", group = "angles", name = "variogram_anis")
plotVariogram <- function(sfe, features, sample_id = "all", color_by = NULL,
                          group = c("none", "sample_id", "features", "angles"),
                          use_lty = TRUE, show_np = TRUE, ncol = NULL,
                          colGeometryName = NULL, annotGeometryName = NULL,
                          reducedDimName = NULL, divergent = FALSE,
                          diverge_center = NULL, swap_rownames = NULL,
                          name = "variogram") {
    rlang::check_installed("gstat")
    group <- match.arg(group)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)

    dfs <- lapply(sample_id, .get_plot_variogram_df, sfe = sfe,
                 features = features, color_by = color_by,
                 colGeometryName = colGeometryName,
                 annotGeometryName = annotGeometryName,
                 reducedDimName = reducedDimName, name = name,
                 show_symbol = !is.null(swap_rownames),
                 swap_rownames = swap_rownames)
    exp_vars <- do.call(rbind, lapply(dfs, function(x) x$exp_vars))
    angles <- sort(unique(exp_vars$dir.hor))
    if (length(angles) > 1L) {
        exp_vars$dir.hor <- factor(as.character(exp_vars$dir.hor),
                                   levels = as.character(angles))
    }
    # set group = "none" when there's only one type in the group
    if (group == "sample_id" && length(sample_id) == 1L ||
        (group == "features" && length(unique(exp_vars$feature)) == 1L) ||
        (group == "angles" && length(angles) == 1L))
        group <- "none"

    var_models <- do.call(rbind, lapply(dfs, function(x) x$var_models))
    model_labels <- do.call(rbind, lapply(dfs, function(x) x$model_labels))
    is_dimred <- is.null(colGeometryName) && is.null(annotGeometryName) &&
        !is.null(reducedDimName) && length(features) > 1L
    if (is_dimred) {
        exp_vars <- .dimred_feature_order(exp_vars)
        var_models <- .dimred_feature_order(var_models)
        model_labels <- .dimred_feature_order(model_labels)
    }
    if (is.data.frame(color_by)) {
        if (length(unique(color_by$cluster)) < 2L) {
            color_by <- NULL
            message("Only one cluster is present. Not using color_by.")
        }
        else {
            names(color_by)[names(color_by) == "cluster"] <- "color_by"
            exp_vars <- merge(exp_vars, color_by, by = c("feature", "sample_id"),
                              all.x = TRUE)
            var_models <- merge(var_models, color_by, by = c("feature", "sample_id"),
                                all.x = TRUE)
        }
    }

    dist <- np <- feature <- dir.hor <- label <- NULL
    base_aes_line <- base_aes_point <- aes(dist, gamma)
    base_aes_np <- aes(dist, gamma, label = np)
    do_color <- !is.null(color_by) | group != "none"
    if (is.null(color_by)) {
        add_aes <- switch (
            group,
            sample_id = aes(color = sample_id),
            features = aes(color = feature),
            angles = aes(color = dir.hor),
            none = list(NULL)
        )
        base_aes_point <- modifyList(base_aes_point, add_aes)
        base_aes_np <- modifyList(base_aes_np, add_aes)
        if (group != "angles") base_aes_line <- modifyList(base_aes_line, add_aes)
        if (group == "angles" || (group == "features" && is_dimred)) {
            name_use <- if (group == "angles") group else if (is_dimred) "component"
            pal <- scale_color_viridis_d(option = "E", end = 0.9,
                                         name = name_use)
        } else {
            pal <- scale_color_manual(values = ditto_colors)
        }
    } else {
        base_aes_line <- base_aes_point <- modifyList(base_aes_line,
                                                      aes(color = color_by))
        base_aes_np <- modifyList(base_aes_np, aes(color = color_by))
        name_use <- if (length(color_by) == 1L && is.atomic(color_by)) color_by
        else if (is.data.frame(color_by)) "cluster"
        else "color_by"
        pal <- .get_pal(exp_vars, list(color = "color_by"), 1,
                        divergent = divergent, diverge_center = diverge_center,
                        name = name_use)
        if (use_lty) {
            add_aes_line <- switch (
                group,
                sample_id = aes(linetype = sample_id),
                features = aes(linetype = feature),
                angles = list(NULL),
                none = list(NULL)
            )
            add_aes_point <- switch (
                group,
                sample_id = aes(shape = sample_id),
                features = aes(shape = feature),
                angles = aes(shape = dir.hor),
                none = list(NULL)
            )
            base_aes_line <- modifyList(base_aes_line, add_aes_line)
            base_aes_point <- modifyList(base_aes_point, add_aes_point)
        } else {
            add_aes_line <- switch (
                group,
                sample_id = aes(group = sample_id),
                features = aes(group = feature),
                angles = list(NULL),
                none = list(NULL)
            )
            base_aes_line <- modifyList(base_aes_line, add_aes_line)
        }
    }
    p <- ggplot() +
        geom_line(data = var_models, mapping = base_aes_line) +
        geom_point(data = exp_vars, mapping = base_aes_point)
    if (show_np)
        p <- p + geom_text(data = exp_vars, mapping = base_aes_np, hjust = 0, vjust = 0)
    if (do_color) p <- p + pal
    # i.e. plotting only one line for one fitted model
    if (group %in% c("none", "angles")) {
        p <- p +
            geom_text(data = model_labels, mapping = aes(label = label),
                      x = max(exp_vars$dist),
                      y = min(min(var_models$gamma), min(exp_vars$gamma)),
                      hjust = 1, vjust = 0)
    }
    v <- c("sample_id", "feature", "dir.hor")
    group2 <- group
    group2[group2 == "angles"] <- "dir.hor"
    group2[group2 == "features"] <- "feature"
    do_facet <- c(length(sample_id) > 1L, length(unique(exp_vars$feature)) > 1L,
                  length(unique(exp_vars$dir.hor)) > 1L)
    # I want angles to be in the columns, consistent with gstat
    facet_by <- sort(setdiff(v[do_facet], group2), decreasing = TRUE)
    if (length(facet_by) == 2L) {
        p <- p +
            facet_grid(rows = vars(!!sym(facet_by[1])),
                       cols = vars(!!sym(facet_by[2])))
    } else if (length(facet_by)) {
        f <- as.formula(paste0("~ ", paste(facet_by, collapse = " + ")))
        p <- p +
            facet_wrap(f, ncol = ncol)
    }
    p <- p +
        labs(x = "Distance", y = "Variance",
             title = "Experimental variogram and fitted variogram model")
    p
}

#' Plot variogram maps
#'
#' Plot variogram maps that show the variogram in all directions in a grid of
#' distances in x and y coordinates.
#'
#' @inheritParams plotVariogram
#' @param plot_np Logical, whether to plot the number of pairs in each distance
#'   bin instead of the variance.
#' @return A ggplot object.
#' @seealso plotVariogram
#' @export
#' @importFrom ggplot2 geom_tile scale_fill_viridis_c
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData()
#' sfe <- colDataUnivariate(sfe, "variogram_map", features = "nCounts",
#' width = 500, cutoff = 5000)
#' plotVariogramMap(sfe, "nCounts")
#'
plotVariogramMap <- function(sfe, features, sample_id = "all", plot_np = FALSE,
                             ncol = NULL,
                             colGeometryName = NULL, annotGeometryName = NULL,
                             reducedDimName = NULL, swap_rownames = NULL,
                             name = "variogram_map") {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    show_symbol <- !is.null(swap_rownames)
    dfs <- lapply(sample_id, function(s) {
        ress <- .get_feature_metadata(
            sfe, features, name, s, colGeometryName,
            annotGeometryName, reducedDimName, show_symbol, swap_rownames
        )
        for (i in seq_along(ress)) {
            ress[[i]]$feature <- names(ress)[[i]]
        }
        out <- do.call(rbind, ress)
        out$sample_id <- s
        out
    })
    df <- do.call(rbind, dfs)
    features <- unique(df$feature)

    dx <- dy <- np.var1 <- var1 <- feature <- NULL
    aes_use <- aes(dx, dy)
    if (plot_np) {
        aes_use <- modifyList(aes_use, aes(fill = np.var1))
        name_use <- "Number\nof pairs"
    } else {
        aes_use <- modifyList(aes_use, aes(fill = var1))
        name_use <- "Variance"
    }
    p <- ggplot(df, aes_use) +
        geom_tile() +
        scale_fill_viridis_c(option = "A", name = name_use) +
        coord_equal()
    facet_inds <- c(length(sample_id) > 1L, length(features) > 1L)
    if (all(facet_inds)) {
        p <- p + facet_grid(rows = vars(sample_id), cols = vars(feature))
    } else if (any(facet_inds)) {
        facet_by <- c("sample_id", "feature")[facet_inds]
        f <- as.formula(paste("~", facet_by))
        p <- p + facet_wrap(f, ncol = ncol)
    }
    p
}

.separate <- function(df, col, delim, names) {
    # Base R implementation of tidyverse separate, to avoid more dependencies
    # Only for the purpose of plotting cross variograms
    ids_split <- strsplit(as.character(df[[col]]), delim)
    ids_split <- lapply(ids_split, function(x) {
        if (length(x) == 1L) rep(x, 2) else x
    })
    ids_split <- do.call(rbind, ids_split)
    colnames(ids_split) <- names
    cbind(df, ids_split)
}

design_ragged_rows <- function(data, rows, cols) {
    # from https://stackoverflow.com/a/75938994/8916916
    facets <- unique(data[c(rows, cols)])
    facet_id <- seq_len(nrow(facets))

    groups <- split(facet_id, facets[rows], drop = TRUE)
    n <- lengths(groups)
    i <- rep(seq_along(n), n)
    j <- sequence(n)

    design <- matrix(NA_integer_, max(i), max(j))
    replace(design, cbind(i, j), facet_id)
}

#' Plot cross variogram
#'
#' Equivalent to \code{gstat::plot.gstatVariogram}, but using ggplot2 to be more
#' customizable.
#'
#' @inheritParams plotVariogram
#' @param res Cross variogram results for one sample, from
#'   \code{\link{calculateBivariate}}. Global bivariate results are not stored
#'   in the SFE object.
#' @return A ggplot object. Unfortunately I haven't figured out a way to collect
#' all the facet labels to the top of the entire plot.
#' @seealso plotCrossVariogramMap
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData()
#' sfe <- sfe[,sfe$in_tissue]
#' sfe <- logNormCounts(sfe)
#'
#' res <- calculateBivariate(sfe, type = "cross_variogram",
#' feature1 = c("Myh1", "Myh2", "Csrp3"), swap_rownames = "symbol")
#' plotCrossVariogram(res)
plotCrossVariogram <- function(res, show_np = TRUE) {
    rlang::check_installed("ggh4x")
    angles <- sort(unique(res$dir.hor))
    if (length(angles) > 1L) {
        res$dir.hor <- factor(as.character(res$dir.hor),
                              levels = as.character(angles))
    }
    res <- .separate(res, "id", delim = "\\.", names = c("id_x", "id_y"))
    l <- names(sort(table(res$id_x), decreasing = TRUE))
    res$id_y <- factor(res$id_y, levels = l)
    res$id_x <- factor(res$id_x, levels = l)
    dist <- gamma <- np <- id_y <- id_x <- dir.hor <- NULL
    aes_use <- aes(dist, gamma)
    if (length(angles) > 1L) aes_use <- modifyList(aes_use, aes(color = dir.hor))
    p <- ggplot(res, aes_use) +
        geom_point() + geom_line()
    if (show_np) p <- p + geom_text(aes(label = np), hjust = 0, vjust = 0)
    if (length(angles) > 1L) p <- p + scale_color_viridis_d(option = "E", end = 0.9,
                                                            name = "Angle")
    design <- design_ragged_rows(res, "id_y", "id_x")
    p <- p + ggh4x::facet_manual(vars(id_y, id_x), design = design,
                                 strip = ggh4x::strip_split(position = c("left", "top")))
    p
}

#' Plot cross variogram map
#'
#' Equivalent to \code{gstat::plot.gstatVariogram}, but using ggplot2 to be more
#' customizable.
#'
#' @inheritParams plotCrossVariogram
#' @inheritParams plotVariogramMap
#' @return A ggplot object.
#' @seealso plotCrossVariogram
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData()
#' sfe <- sfe[,sfe$in_tissue]
#' sfe <- logNormCounts(sfe)
#'
#' res <- calculateBivariate(sfe, type = "cross_variogram_map",
#' feature1 = c("Myh1", "Myh2", "Csrp3"), swap_rownames = "symbol",
#' width = 500, cutoff = 2000)
#' plotCrossVariogramMap(res)
#'
plotCrossVariogramMap <- function(res, plot_np = FALSE) {
    rlang::check_installed("ggh4x")
    cols_use <- setdiff(names(res), c("dx", "dy"))
    is_np <- grepl("^np", x = cols_use)
    if (plot_np) {
        cols_use <- cols_use[is_np]
        res <- res[,c(cols_use, "dx", "dy")]
        names(res) <- gsub("^np\\.", "", names(res))
        name_use <- "Number\nof pairs"
    } else {
        cols_use <- cols_use[!is_np]
        res <- res[,c(cols_use, "dx", "dy")]
        name_use <- "Variance"
    }
    res <- reshape(res, varying = cols_use, direction = "long",
                   v.names = "value", timevar = "variable",
                   times = cols_use)
    res <- .separate(res, "variable", delim = "\\.", names = c("id_x", "id_y"))
    l <- names(sort(table(res$id_x), decreasing = TRUE))
    res$id_y <- factor(res$id_y, levels = l)
    res$id_x <- factor(res$id_x, levels = l)
    design <- design_ragged_rows(res, "id_y", "id_x")
    dx <- dy <- value <- id_y <- id_x <- NULL
    ggplot(res, aes(dx, dy, fill = value)) +
        geom_tile() + scale_fill_viridis_c(option = "A", name = name_use) +
        ggh4x::facet_manual(vars(id_y, id_x), design = design,
                            strip = ggh4x::strip_split(position = c("left", "top")))
}
