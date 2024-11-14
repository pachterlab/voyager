# 10. What to do with the image when using geom_sf
# To do:
# 2. reverse_y? It's kind of hard to do that with sf.

#' Get beginning and end of palette to center a divergent palette
#'
#' This function is no longer used internally as it's unnecessary for
#' \code{scico} divergent palettes. But it can be useful when using divergent
#' palettes outside \code{scico} where one must specify beginning and end but
#' not midpoint, to override the default palette.
#'
#' @param values Numeric vector to be colored.
#' @param diverge_center Value to center on, defaults to 0.
#' @return A numeric vector of length 2, the first element is for beginning, and
#' the second for end. The values are between 0 and 1.
#' @concept Spatial plotting
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
        shape <- fixed[["shape"]]
        if (!is.null(shape) && shape > 20L) {
            names_use <- c(names_use, "fill")
        }
    } else if (type %in% c("LINESTRING", "MULTILINESTRING")) {
        if (isTRUE(all.equal(0, fixed$linewidth))) fixed$linewidth <- 0.5
        names_use <- c("linewidth", "linetype", "alpha", "color")
    } else {
        # i.e. polygons
        names_use <- c("linewidth", "linetype", "fill", "color", "alpha")
        if (isTRUE(all.equal(0, fixed$linewidth))) fixed$color <- NA
    }
    fixed_applicable <- .drop_null_list(fixed[names_use])
    fixed_applicable
}

.get_pal <- function(df, feature_aes, option, divergent, diverge_center,
                     name = waiver(), dark = FALSE) {
    feature_aes <- feature_aes[names(feature_aes) %in% c("fill", "color", "z")]
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
        pal <- pal_fun(values = .pal, na.value = "gray", name = name)
    } else {
        if (divergent) {
            if (dark) {
                .pal <- switch(option,
                               "berlin",
                               "vanimo"
                )
            } else {
                .pal <- switch(option,
                               "roma",
                               "bam"
                )
            }
            pal_fun <- switch(.aes,
                fill = scale_fill_scico,
                color = scale_color_scico,
                z = scale_fill_scico
            )
            .direction <- if (.pal == "berlin") 1 else -1
            pal <- pal_fun(
                palette = .pal, direction = .direction, midpoint = diverge_center,
                na.value = "gray", name = name
            )
        } else {
            if (dark) {
                .pal <- switch(option,
                               "nuuk",
                               "acton")
                pal_fun <- switch(.aes,
                                   fill = scale_fill_scico,
                                   color = scale_color_scico,
                                   z = scale_fill_scico
                )
                pal <- pal_fun(na.value = "gray", name = name, palette = .pal)
            } else {
                pal_fun <- switch(.aes,
                                  fill = scale_fill_distiller,
                                  color = scale_color_distiller,
                                  z = scale_fill_distiller
                )
                .pal <- switch(option,
                               "Blues",
                               "YlOrRd"
                )
                pal <- pal_fun(na.value = "gray", palette = .pal, direction = 1,
                               name = name)
            }
        }
    }
    pal
}

.fill_defaults <- function(fixed) {
    defaults <- list(
        size = 0, linewidth = 0, shape = 16,
        linetype = 1, alpha = 1, color = "black", fill = "gray80",
        divergent = FALSE, diverge_center = NULL
    )
    fill <- defaults[setdiff(names(defaults), names(fixed))]
    out <- .drop_null_list(c(fixed, fill))
    if (is.na(out$fill) && "linewidth" %in% names(fill)) {
        out$linewidth <- 0.5
    }
    out
}

#' @importFrom ggplot2 element_rect element_text theme margin %+replace% element_line
#' theme_gray rel
.dark_theme <- function(show_axes = FALSE) {
    # From Seurat but no axes
    black.background <- element_rect(fill = 'black')
    black.background.no.border <- element_rect(fill = 'black', linewidth = 0)
    font.margin <- 4
    white.text <- element_text(
        colour = 'white',
        margin = margin(
            t = font.margin,
            r = font.margin,
            b = font.margin,
            l = font.margin
        )
    )
    # Create the dark theme
    if (show_axes) {
        theme_gray() %+replace%
            theme(axis.text = element_text(size = rel(0.8), colour = "white"),
                  panel.border = element_rect(fill = NA, colour = "grey80"),
                  panel.grid = element_line(colour = "gray30"),
                  axis.ticks = element_line(colour = "grey30"),
                  panel.grid.minor = element_line(linewidth = rel(0.5)),
                  plot.background = black.background,
                  panel.background = black.background,
                  legend.background = black.background,
                  legend.box.background = black.background.no.border,
                  legend.key = black.background.no.border,
                  strip.background = black.background,
                  text = white.text,
                  validate = TRUE)
    } else {
        theme_void() %+replace%
            theme(
                #   Set background colors
                plot.background = black.background,
                panel.background = black.background,
                legend.background = black.background,
                legend.box.background = black.background.no.border,
                legend.key = black.background.no.border,
                strip.background = black.background,
                #   Set text colors
                text = white.text,
                #   Validate the theme
                validate = TRUE
            )
    }
}

#' @importFrom sf st_drop_geometry st_geometry_type
#' @importFrom ggplot2 ggplot geom_sf scale_fill_manual
#' scale_color_manual scale_fill_distiller scale_color_distiller geom_polygon
#' geom_segment stat_density2d waiver stat_summary_2d stat_summary_hex
#' @importFrom scico scale_fill_scico scale_color_scico
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom rlang syms !!!
.plot_var_sf <- function(df, annot_df, img_df, channel, type, type_annot, feature_aes,
                         feature_fixed, annot_aes, annot_fixed, tx_fixed, divergent,
                         diverge_center,annot_divergent, annot_diverge_center,
                         ncol_sample, scattermore, pointsize,
                         bins, summary_fun, hex, maxcell, show_axes, dark, palette,
                         tx_df) {
    # Add annotGeometry if present
    if (!is.null(annot_df)) {
        annot_fixed <- .get_applicable(type_annot, annot_fixed)
        annot_fixed <- annot_fixed[setdiff(names(annot_fixed), names(annot_aes))]
        if ("color" %in% names(annot_aes) && annot_fixed$linewidth == 0) {
            annot_fixed$linewidth <- 0.5
        }
        geom_annot <- do.call(geom_sf, c(
            list(mapping = aes(!!!syms(annot_aes)), data = annot_df),
            annot_fixed
        ))
        pal_annot <- .get_pal(annot_df, annot_aes, 2, annot_divergent,
                              annot_diverge_center, dark = dark)
    }
    # Deal with gene symbols that are not legal R object names
    df_names_orig <- names(df)
    names(df) <- make.names(names(df))
    feature_aes <- lapply(feature_aes, make.names)
    inds <- !df_names_orig %in% names(df)
    ind_plot <- which(names(df) == unlist(feature_aes)) # should be only one
    if (any(inds) && length(ind_plot))
        name_show <- df_names_orig[ind_plot] else name_show <- waiver()

    p <- ggplot()
    data <- NULL
    if (!is.null(img_df)) {
        # Check if it's RGB
        img <- img_df$data[[1]]
        img_dim <- dim(img)
        p <- p + geom_spi_rgb(data = img_df, aes(spi = data),
                              palette = palette)
    }
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
    if (scattermore) {
        geom_use <- do.call(scattermore::geom_scattermore,
                            c(list(mapping = aes(!!!syms(c(list(x = "X", y = "Y"),
                                                           feature_aes))),
                                   data = df, pointsize = pointsize),
                              feature_fixed))
    } else if (!is.null(bins)) {
        names(feature_aes)[names(feature_aes) == "color"] <- "z"
        name_show <- paste("Aggregated", feature_aes[["z"]], sep = "\n")
        feature_fixed <- feature_fixed["alpha"]
        stat_fun <- if (hex) stat_summary_hex else stat_summary_2d
        geom_use <- do.call(stat_fun,
                            c(list(mapping = aes(!!!syms(c(list(x = "X", y = "Y"),
                                                           feature_aes))),
                                   data = df, bins = bins, fun = summary_fun),
                              feature_fixed))
    } else {
        geom_use <- do.call(geom_sf, c(
            list(mapping = aes(!!!syms(feature_aes)), data = df),
            feature_fixed
        ))
    }
    p <- p + geom_use
    if (scattermore || !is.null(bins)) p <- p + coord_equal()
    # Palette
    pal <- .get_pal(df, feature_aes, 1, divergent, diverge_center,
                    name = name_show, dark = dark)
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
    if (!is.null(tx_df)) {
        tx_fixed <- .get_applicable("MULTIPOINT", tx_fixed)
        if (length(unique(tx_df$gene)) > 1L) {
            tx_fixed <- tx_fixed[names(tx_fixed) != "shape"]
            geom_tx <- do.call(geom_sf, c(
                list(mapping = aes(shape = gene), data = tx_df),
                tx_fixed
            ))
        } else geom_tx <- do.call(geom_sf, c(
            list(data = tx_df), tx_fixed
        ))
            p <- p + geom_tx
    }
    if ("sample_id" %in% names(df) && length(unique(df$sample_id)) > 1L) {
        p <- p +
            facet_wrap(~sample_id, ncol = ncol_sample)
    }

    if (dark) {
        p <- p + .dark_theme(show_axes)
    } else {
        if (show_axes) p <- p + theme_bw() else p <- p + theme_void()
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

.wrap_spatial_plots <- function(df, annot_df, img_df, channel, type_annot, values, aes_use,
                                annot_aes, annot_fixed, tx_fixed, size, shape, linewidth,
                                linetype, alpha, color, fill, ncol, ncol_sample,
                                divergent, diverge_center, annot_divergent,
                                annot_diverge_center, scattermore, pointsize,
                                bins, summary_fun, hex, maxcell, show_axes, dark,
                                palette, tx_df, rowGeometryFeatures, ...) {
    feature_fixed <- list(
        size = size, linewidth = linewidth, shape = shape, linetype = linetype,
        alpha = alpha, color = color, fill = fill
    )
    type <- if(scattermore || !is.null(bins)) "POINT" else .get_generalized_geometry_type(df)
    features <- names(values)
    if (!is.null(rowGeometryFeatures) && all(rowGeometryFeatures %in% features)) {
        tx_dfs <- split(tx_df, tx_df$gene)
        plots <- lapply(features, function(n) {
            feature_aes_name <- .get_feature_aes(df[[n]], type, aes_use, shape)
            feature_aes <- setNames(list(n), feature_aes_name)
            .plot_var_sf(
                df, annot_df, img_df, channel, type, type_annot, feature_aes, feature_fixed,
                annot_aes, annot_fixed, tx_fixed, divergent, diverge_center,
                annot_divergent, annot_diverge_center, ncol_sample, scattermore,
                pointsize, bins, summary_fun, hex, maxcell, show_axes, dark, palette,
                tx_dfs[[n]]
            )
        })
    } else {
        plots <- lapply(names(values), function(n) {
            feature_aes_name <- .get_feature_aes(df[[n]], type, aes_use, shape)
            feature_aes <- setNames(list(n), feature_aes_name)
            .plot_var_sf(
                df, annot_df, img_df, channel, type, type_annot, feature_aes, feature_fixed,
                annot_aes, annot_fixed, tx_fixed, divergent, diverge_center,
                annot_divergent, annot_diverge_center, ncol_sample, scattermore,
                pointsize, bins, summary_fun, hex, maxcell, show_axes, dark, palette,
                tx_df
            )
        })
    }
    if (length(plots) > 1L) {
        out <- wrap_plots(plots, ncol = ncol, ...)
        if (dark) {
            out <- out +
                plot_annotation(theme = theme(text = element_text(color = "white"),
                                              plot.background = element_rect(fill = "black")))
        }
    } else {
        out <- plots[[1]]
    }
    out
}

#' @importFrom sf st_as_sfc st_bbox st_intersection
.bbox_sample <- function(df, bbox) {
    if (is(df, "sf")) {
        # Only for one sample
        bbox_use <- st_as_sfc(st_bbox(bbox))
        suppressWarnings(df <- st_intersection(df, bbox_use))
        #df <- df[bbox_use,]
        df$geometry <- st_geometry(df) - bbox[c("xmin", "ymin")]
    } else {
        # i.e. centroids. Assume column names x and y
        inds <- df$x > bbox["xmin"] & df$x < bbox["xmax"] &
            df$y > bbox["ymin"] & df$y < bbox["ymax"]
        df <- df[inds,, drop = FALSE]
        # Also need to subtract xmin and ymin
        df$x <- df$x - bbox["xmin"]
        df$y <- df$y - bbox["ymin"]
    }
    if (nrow(df) == 0L)
        stop("The bounding box does not overlap with the geometry.")
    df
}

# Burning question: what if the bboxes of different samples are very different
# so there will be a lot of empty space in the facetted plot?
# I suppose I'll subtract xmin and ymin. The original coordinates are known
# from the bbox anyway.
# But what if bboxes of different samples have very different sizes?
# Well, it's up to the user to not to do it.
.crop <- function(df, bbox) {
    if (is.null(bbox)) return(df)
    if (!is.atomic(bbox) && !is.matrix(bbox)) {
        stop("bbox must be either a numeric vector or a matrix.")
    }
    names_use <- c("xmin", "ymin", "xmax", "ymax")
    if (is.matrix(bbox)) {
        if (nrow(bbox) == 1L) bbox <- setNames(bbox, colnames(bbox))
        if (ncol(bbox) == 1L) bbox <- setNames(bbox, rownames(bbox))
    }
    if (is.vector(bbox)) {
        if (!is.numeric(bbox)) stop("bbox must be a numeric vector")
        if (length(bbox) != 4L)
            stop("bbox must have length 4, corresponding to xmin, ymin, xmax, and ymax.")

        if (!is.null(names(bbox)) && !setequal(names(bbox), names_use)) {
            stop("Names of bbox must be same as the set ",
                 paste(names_use, collapse = ", "))
        }
        if (is.null(names(bbox))) {
            warning("No names available for bbox. Assuming ",
                    paste(names_use, collapse = ", "), ", in this order.")
            names(bbox) <- names_use
        }
    }
    if (is.matrix(bbox)) {
        if (setequal(colnames(bbox), names_use)) bbox <- t(bbox)
        if (nrow(bbox) != 4L)
            stop("bbox must have 4 rows, corresponding to xmin, ymin, xmax, and ymax.")
        if (!setequal(rownames(bbox), names_use))
            stop("Row names of bbox must be same as the set ",
                 paste(names_use, collapse = ", "))
        if (length(setdiff(unique(df$sample_id), colnames(bbox)))) {
            stop("Column names of bbox must match the sample IDs")
        }
        if (!"sample_id" %in% names(df)) {
            warning("Only the first column of the matrix bbox will be used.")
            bbox <- bbox[,1]
        } else
            bbox <- bbox[,unique(df$sample_id)]
    }

    if (!"sample_id" %in% names(df) || length(unique(df$sample_id)) == 1L) {
        df <- .bbox_sample(df, bbox)
    } else {
        df_split <- split(df, df$sample_id)
        if (is.vector(bbox)) {
            df_split <- lapply(df_split, .bbox_sample, bbox = bbox)
        } else {
            samples <- names(df_split)
            df_split <- lapply(samples, function(n)
                .bbox_sample(df_split[[n]], bbox[,n]))
        }
        df <- do.call(rbind, df_split)
    }
    df
}

#' @importFrom SpatialFeatureExperiment isFull imgSource getPixelSize ext
#' aggBboxes cropImg imageIDs translateImg toSpatRasterImage getParams
#' @importFrom sf st_area
#' @importFrom terra RGB<- rast
.find_res <- function(bfi, maxcell) {
    check_installed("RBioFormats")
    coreMetadata <- RBioFormats::coreMetadata
    metas <- RBioFormats::read.metadata(imgSource(bfi))
    n_series <- RBioFormats::seriesCount(metas)
    cms <- coreMetadata(metas)
    if (n_series == 1L) return(1L)

    if (isFull(bfi)) {
        dims <- data.frame(x = vapply(cms, function(x) x$sizeX, FUN.VALUE = numeric(1)),
                           y = vapply(cms, function(x) x$sizeY, FUN.VALUE = numeric(1)))
        ncells <- dims$x*dims$y
    } else {
        # Get the number of pixels within the extent at each resolution
        ncells <- vapply(seq_len(n_series), function(i) {
            ps <- getPixelSize(imgSource(bfi), resolution = i)
            psx <- ps[1]; psy <- ps[2]
            bb <- ext(bfi)
            npx_x <- (bb["xmax"] - bb["xmin"])/psx
            npx_y <- (bb["ymax"] - bb["ymin"])/psy
            npx_x*npx_y
        }, FUN.VALUE = numeric(1))
    }
    n_use <- min(which(ncells < 1.1*maxcell)) # 1.1 as in resample_spat
    return(n_use)
}

.get_n_channels <- function(img) {
    d <- dim(img)
    if (length(d) == 2L) return(1L)
    as.integer(d[[3]])
}
# What I want: option to assign individual grayscale images to different channels
# Option to select one or more channels from a multi-channel image
.combine_channels <- function(imgs, channel_assign) {
    # Already checked, cropped, and converted to SpatRaster
    # BioFormatsImage
    # choose highest resolution with fewer than maxcell when multiple res are present
    # Convert to spi
    # ExtImage: convert to spi so can use custom color map in terra::as.raster
    # Combine different image classes: all use spi
    # For testing, need to use xenium v2 example data, everything is connected
    # All convert to spi, find the one with the lowest resolution, then
    # resample all the others to the same resolution before combining them into
    # rgb channels
    ncells <- vapply(imgs, ncell, FUN.VALUE = numeric(1))
    dims <- vapply(imgs, function(img) dim(img)[1:2], FUN.VALUE = numeric(2))
    same_dims <- length(unique(dims[1,])) == 1L & length(unique(dims[2,])) == 1L
    # TODO: What if different channels have different extents
    if (!same_dims) {
        ind <- which.min(ncells)
        ind_resample <- seq_along(imgs)[-ind]
        imgs[ind_resample] <- lapply(imgs[ind_resample], function(img) {
            # It's just for plotting so I don't care about the method used to resample
            terra::resample(img, imgs[[ind]])
        })
    }
    ch <- c("r", "g", "b")
    if (length(imgs) == 3L) {
        names(imgs) <- channel_assign
        imgs <- imgs[ch]
        out <- rast(imgs)
    } else if (length(imgs) == 2L) {
        # Create raster of all 0's. Take names from the channels
        bl <- rast(imgs[[1]], vals = 0)
        bl_channel <- setdiff(ch, channel_assign)
        ol <- setNames(c(imgs, bl), c(channel_assign, bl_channel))
        ol <- ol[ch]
        out <- rast(ol)
    }
    RGB(out) <- 1:3
    return(out)
}

.subset_channels <- function(img, channel) {
    # Mainly to assign the RGB channels
    ch <- c("r", "g", "b")
    if (length(channel) == 1L) {
        return(img[[channel]])
    } else if (length(channel) == 2L) {
        bl <- rast(img, nlyrs = 1, vals = 0)
        bl_channel <- setdiff(ch, channel)
        chs <- c(names(channel), bl_channel)
        ch_ord <- match(ch, chs)
        out <- c(img[[channel]], bl)
        out <- out[[ch_ord]]
    } else if (length(channel) == 3L) {
        if (!is.null(names(channel))) {
            channel <- channel[ch]
        }
        out <- img[[channel]]
    }
    RGB(out) <- 1:3
    return(out)
}

.normalize_channels <- function(img) {
    max_values <- terra::minmax(img, compute = TRUE)[2,]
    for (i in seq_len(nlyr(img))) {
        if (max_values[[i]] > 0) img[[i]] <- img[[i]]/max_values[[i]]
    }
    img
}

#' @importFrom memuse Sys.meminfo
.get_img_df_sample <- function(sample_id, df, image_id, channel, bbox, maxcell,
                               normalize_channels) {
    # For each sample
    df <- df[df$sample_id == sample_id,]
    image_id <- image_id[match(df$image_id, image_id)]
    if (!is.null(channel)) {
        if (!is.numeric(channel) || any(channel < 1L))
            stop("channel must be numeric indices")
        if (length(channel) > 3L) {
            stop("Only up to 3 channels can plotted at once in an RGB image")
        }
        if (!is.null(names(channel)) && !any(names(channel) %in% c("r", "g", "b"))) {
            stop("Names of channel indices must be among 'r', 'g', 'b'")
        }
        if (length(channel) == 2L && is.null(names(channel))) {
            stop("channel of length 2 must have names to specify which of the RGB channels to use")
        }
        if (length(image_id) > 1L) {
            warning("Cannot use multiple images as different channels when ",
                    "argument channel is specified to select channels in one image. ",
                    "Only using the first image.")
            df <- df[1,,drop = FALSE]
            image_id <- image_id[1]
        }
        n_channels <- .get_n_channels(df$data[[1]])
        if (any(channel > n_channels))
            stop("channel index out of bound")
    }
    imgs <- df$data
    n_channels <- vapply(imgs, .get_n_channels, FUN.VALUE = integer(1L))
    if (length(image_id) > 1L) {
        # Check names if length < 3L, don't use red + green by default
        if (length(image_id) == 2L && is.null(names(image_id))) {
            stop("image_id of length 2 must have names to specify which of the RGB channels to use")
        }
        if (length(imgs) > 3L)
            stop("Colorization allows up to 3 channels")
        # All images must have only 1 channel
        if (any(n_channels > 1L)) {
            ind <- which(n_channels > 1L)
            stop("All images to be combined as different channels must only have 1 channel. ",
                 "Image(s) number ", paste(ind, collapse = ", "), " has/have ",
                 "multiple channels.")
        }
    } else if (n_channels > 1L && is.null(channel)) {
        if (n_channels == 3L) channel <- 1:3
        else stop("Argument `channel` must be specified to map channels to colors.")
    }
    # Crop
    if (!is.null(bbox)) {
        # For SpatRaster, for huge images on disk, downsample before cropping if
        # the bbox is a large part of the image that also needs to be written to
        # disk
        imgs <- lapply(imgs, function(x) {
            if (is(x, "SpatRasterImage")) {
                tot_area <- ext(x) |> st_bbox() |> st_as_sfc() |> st_area()
                bb_area <- bbox |> st_bbox() |> st_as_sfc() |> st_area()
                bb_prop <- bb_area/tot_area
                if (!terra::inMemory(x)) {
                    # Shouldn't need to write the cropped part to disk just for a plot
                    tot_size <- file.info(imgSource(x))$size
                    mem_free <- Sys.meminfo()$freeram |> as.numeric()
                    if (bb_prop * tot_size > mem_free/2) {
                        # The /2 since images take more RAM than disk space when compressed
                        # But what if maxcell_tot is still too large?
                        # What if it's larger than the fullres cropped area?
                        maxcell_tot <- maxcell/bb_prop
                        ds_prop <- maxcell_tot/ncell(x)
                        if (ds_prop < bb_prop)
                            x@image <- resample_spat(x@image, maxcell_tot)
                    }
                }
            }
            out <- cropImg(x, bbox)
            translateImg(out, -bbox[c("xmin", "ymin")])
        })
    }
    # All convert to SpatRaster
    imgs <- lapply(imgs, function(img) {
        if (is(img, "BioFormatsImage")) {
            res_use <- .find_res(img, maxcell)
            spi <- toSpatRasterImage(img, resolution = res_use, save_geotiff = FALSE)
        } else if (is(img, "ExtImage")) {
            spi <- toSpatRasterImage(img, save_geotiff = FALSE)
        } else spi <- img
        spi |> resample_spat(maxcell)
    })
    # Combine channels
    if (length(image_id) > 1L)
        img <- .combine_channels(imgs, names(image_id))
    # Subset channels
    else if (!is.null(channel)) {
        img <- .subset_channels(imgs[[1]], channel)
    } else img <- imgs[[1]]
    if (normalize_channels) img <- .normalize_channels(img)
    # Output: should have only 1 image left
    return(img)
}

.get_img_df <- function(sfe, sample_id, image_id, channel, bbox, maxcell,
                        normalize_channels) {
    image_id <- image_id %||% imageIDs(sfe)[1]
    img_df <- imgData(sfe)
    img_df <- img_df[img_df$sample_id %in% sample_id & img_df$image_id %in% image_id,
                     c("sample_id", "data", "image_id")]

    # Edge case: when different images are used for different channels and there're
    # multiple samples, some samples don't have images for all the channels
    imgs <- lapply(sample_id, .get_img_df_sample, df = img_df,
                   image_id = image_id, channel = channel, bbox = bbox,
                   maxcell = maxcell, normalize_channels = normalize_channels)
    data.frame(sample_id = sample_id, data = I(imgs))
}

#' @importFrom rlang check_installed
.plotSpatialFeature <- function(sfe, values, colGeometryName, sample_id, ncol,
                                ncol_sample, annotGeometryName, annot_aes,
                                annot_fixed, tx_fixed, bbox, image_id, channel, aes_use, divergent,
                                diverge_center, annot_divergent,
                                annot_diverge_center, size, shape, linewidth,
                                linetype, alpha, color, fill, scattermore,
                                pointsize, bins, summary_fun, hex, maxcell,
                                show_axes, dark, palette, normalize_channels,
                                rowGeometryName, rowGeometryFeatures, tx_file,...) {
    df <- colGeometry(sfe, colGeometryName, sample_id = sample_id)
    df$sample_id <- colData(sfe)$sample_id[colData(sfe)$sample_id %in% sample_id]
    # In case of illegal names
    names_orig <- names(values)
    df <- cbind(df[,c("geometry", "sample_id")], values)
    df <- .crop(df, bbox)
    names(df)[!names(df) %in% c("geometry", "sample_id")] <- names_orig
    type_df <- .get_generalized_geometry_type(df)
    if (type_df %in% c("POLYGON", "MULTIPOLYGON") && is.na(fill) && linewidth == 0) {
        linewidth <- 0.3
    }
    if (scattermore || !is.null(bins)) {
        if (scattermore)
            check_installed("scattermore",
                            reason = "to plot points with scattermore")
        if (type_df != "POINT") {
            message("scattermore and binning only apply to points. Using centroids.")
            df_coords <- as.data.frame(st_coordinates(st_centroid(st_geometry(df))))
        } else {
            df_coords <- as.data.frame(st_coordinates(df))
        }
        df <- cbind(st_drop_geometry(df), df_coords)
    }
    # Will use separate ggplots for each feature so each can have its own color scale
    if (!is.null(annotGeometryName)) {
        annot_df <- annotGeometry(sfe, annotGeometryName, sample_id)
        type_annot <- .get_generalized_geometry_type(annot_df)
        annot_df <- .crop(annot_df[,c(unlist(annot_aes), "sample_id")], bbox)
    } else {
        annot_df <- NULL
        type_annot <- NULL
    }
    if (!is.null(rowGeometryName)) {
        tx_df <- .get_tx_df(sfe, data_dir = NULL, tech = NULL, file = NULL,
                            sample_id = sample_id, spatialCoordsNames = c("X", "Y"),
                            gene_col = "gene", bbox = bbox, gene = rowGeometryFeatures,
                            return_sf = TRUE, rowGeometryName = rowGeometryName,
                            geoparquet_file = tx_file)
    } else tx_df <- NULL
    if (!is.null(image_id)) {
        img_df <- .get_img_df(sfe, sample_id, image_id, channel, bbox, maxcell,
                              normalize_channels)
    } else img_df <- NULL
    if (is(img_df, "DataFrame") && !nrow(img_df)) img_df <- NULL
    .wrap_spatial_plots(
        df, annot_df, img_df, channel, type_annot, values, aes_use,
        annot_aes, annot_fixed, tx_fixed, size, shape, linewidth, linetype, alpha,
        color, fill, ncol, ncol_sample, divergent,
        diverge_center, annot_divergent, annot_diverge_center, scattermore,
        pointsize, bins, summary_fun, hex, maxcell, show_axes, dark, palette,
        tx_df, rowGeometryFeatures,...
    )
}

.getRowGeometryFeatures <- function(sfe, features, rowGeometryName, rowGeometryFeatures) {
    if (!is.null(rowGeometryName)) {
        if (is.null(rowGeometryFeatures)) {
            rowGeometryFeatures <- intersect(features, rownames(sfe))
        } else {
            rowGeometryFeatures <- intersect(rowGeometryFeatures, rownames(sfe))
        }
        if (!length(rowGeometryFeatures)) rowGeometryName <- NULL
    }
    list(name = rowGeometryName, features = rowGeometryFeatures)
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
#' In the light theme, for continuous variables, the Blues palette from
#' colorbrewer is used if \code{divergent = FALSE}, and the roma palette from
#' the \code{scico} package if \code{divergent = TRUE}. In the dark theme, the
#' nuuk palette from \code{scico} is used if \code{divergent = FALSE}, and the
#' berlin palette from \code{scico} is used if \code{divergent = TRUE}. For
#' discrete variables, the \code{dittoSeq} palette is used.
#'
#' For annotation, the YlOrRd colorbrewer palette is used for continuous
#' variables in the light theme. In the dark theme, the acton palette from
#' \code{scico} is used when \code{divergent = FALSE} and the vanimo palette
#' from \code{scico} is used when \code{divergent = FALSE}. The other end of the
#' \code{dittoSeq} palette is used for discrete variables.
#'
#' Each individual palette should be colorblind friendly, but when plotting
#' continuous variables coloring a \code{colGeometry} and a \code{annotGeometry}
#' simultaneously, the combination of the two palettes is not guaranteed to be
#' colorblind friendly.
#'
#' In addition, when plotting an image behind the geometries, the colors of the
#' image may distort color perception of the values of the geometries.
#'
#' \code{theme_void} is used for all spatial plots in this package, because the
#' units in the spatial coordinates are often arbitrary. This can be overriden
#' to show the axes by using a different theme as normally done in
#' \code{ggplot2}.
#'
#' @inheritParams plotCorrelogram
#' @inheritParams calculateUnivariate
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param features Features to plot, must be in rownames of the gene count
#'   matrix, colnames of colData or a colGeometry.
#' @param bbox A bounding box to specify a smaller region to plot, useful when
#'   the dataset is large. Can be a named numeric vector with names "xmin",
#'   "xmax", "ymin", and "ymax", in any order. If plotting multiple samples, it
#'   should be a matrix with sample IDs as column names and "xmin", "ymin",
#'   "xmax", and "ymax" as row names. If multiple samples are plotted but
#'   \code{bbox} is a vector rather than a matrix, then the same bounding box
#'   will be used for all samples. You may see points at the edge of the
#'   geometries if the intersection between the bounding box and a geometry
#'   happens to be a point there. If \code{NULL}, then the entire tissue is
#'   plotted.
#' @param divergent Logical, whether a divergent palette should be used.
#' @param diverge_center If \code{divergent = TRUE}, the center from which the
#'   palette should diverge. If \code{NULL}, then not centering.
#' @param size Fixed size of points. For points defaults to 0.5. Ignored if
#'   \code{size_by} is specified.
#' @param shape Fixed shape of points, ignored if \code{shape_by} is specified
#'   and applicable.
#' @param linewidth Width of lines, including outlines of polygons. For
#'   polygons, this defaults to 0, meaning no outlines.
#' @param linetype Fixed line type, ignored if \code{linetype_by} is specified
#'   and applicable.
#' @param color Fixed color for \code{colGeometry} if \code{color_by} is not
#'   specified or not applicable, or for \code{annotGeometry} if
#'   \code{annot_color_by} is not specified or not applicable.
#' @param fill Similar to \code{color}, but for fill.
#' @param alpha Transparency.
#' @param exprs_values Integer scalar or string indicating which assay of x
#'   contains the expression values.
#' @param annotGeometryName Name of a \code{annotGeometry} of the SFE object, to
#'   annotate the gene expression plot.
#' @param rowGeometryName Name of a \code{rowGeometry} of the SFE object to
#'   plot.
#' @param rowGeometryFeatures Which features from \code{rowGeometry} to plot.
#'   Can only be a small number to avoid overplotting. Different features are
#'   distinguished by point shape. By default (\code{NULL}), when
#'   \code{rowGeometryName} is specified, this will be whichever items in
#'   \code{features} that are also in the row names of the SFE object. If
#'   features specified for this argument are not the same as or a subset of
#'   those in argument \code{features}, then the spots of all features specified
#'   here will be plotted, differentiated by point shape.
#' @param annot_aes A named list of plotting parameters for the annotation sf
#'   data frame. The names are which geom (as in ggplot2, such as color and
#'   fill), and the values are column names in the annotation sf data frame.
#'   Tidyeval is NOT supported.
#' @param annot_fixed Similar to \code{annot_aes}, but for fixed aesthetic
#'   settings, such as \code{color = "gray"}. The defaults are the same as the
#'   relevant defaults for this function.
#' @param tx_fixed Similar to \code{annot_fixed}, but to specify fixed aesthetic
#'   for transcript spots.
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
#' @param scattermore Logical, whether to use the \code{scattermore} package to
#'   greatly speed up plotting numerous points. Only used for POINT
#'   \code{colGeometries}. If the geometry is not POINT, then the centroids are
#'   used. Recommended for plotting hundreds of thousands or more cells where
#'   the cell polygons can't be seen when plotted due to the large number of
#'   cells and small plot size such as when plotting multiple panels for
#'   multiple features.
#' @param bins If binning the \code{colGeometry} in space due to large number of
#'   cells or spots, the number of bins, passed to \code{\link{geom_bin2d}} or
#'   \code{\link{geom_hex}}. If \code{NULL} (default), then the
#'   \code{colGeometry} is plotted without binning. If binning, a point geometry
#'   is recommended. If the geometry is not point, then the centroids will be
#'   used.
#' @param summary_fun Function to summarize the feature value when the
#'   \code{colGeometry} is binned.
#' @param hex Logical, whether to use \code{\link{geom_hex}}. Note that
#'   \code{geom_hex} is broken in \code{ggplot2} version 3.4.0. Please update
#'   \code{ggplot2} if you are getting horizontal stripes when \code{hex =
#'   TRUE}.
#' @param image_id ID of the image to plot behind the geometries. If
#'   \code{NULL}, then not plotting images. Use \code{\link{imgData}} to see
#'   image IDs present. To plot multiple grayscale images as different RGB
#'   channels, use a named vector here, whose names are channel names (r, g, b),
#'   and values are image_ids of the corresponding images. The RGB colorization
#'   may not be colorblind friendly. When plotting multiple samples, it is
#'   assumed that the same image_id is used for each channel across different
#'   samples.
#' @param channel Numeric vector indicating which channels in a multi-channel
#'   image to plot. If \code{NULL}, grayscale is plotted if there is 1 channel
#'   and RGB for the first 3 channels. The numeric vector can be named (r, g, b)
#'   to indicate which channel maps to which color. The RGB colorization may not
#'   be colorblind friendly. This argument cannot be specified while
#'   \code{image_id} is a named vector to plot different grayscale images as
#'   different channels.
#' @param maxcell Maximum number of pixels to plot in the image. If the image is
#'   larger, it will be resampled so it have less than this number of pixels to
#'   save memory and for faster plotting. We recommend reducing this number when
#'   plotting multiple facets.
#' @param pointsize Radius of rasterized point in \code{scattermore}. Default to
#'   0 for single pixels (fastest).
#' @param dark Logical, whether to use dark theme. When using dark theme, the
#'   palette will have lighter color represent higher values as if glowing in
#'   the dark. This is intended for plotting gene expression on top of
#'   fluorescent images.
#' @param normalize_channels Logical, whether to normalize each channel of the
#'   image individually. Should be \code{FALSE} for bright field color images
#'   such as H&E but should set to \code{TRUE} for fluorescent images.
#' @param tx_file File path to GeoParquet file of the transcript spots if you
#'   don't wish to load all transcript spots into the SFE object. See
#'   \code{\link{formatTxSpots}} on generating such a GeoParquet file.
#' @param palette Vector of colors to use to color grayscale images.
#' @param show_axes Logical, whether to show axes.
#' @param ... Other arguments passed to \code{\link{wrap_plots}}.
#' @return A \code{ggplot2} object if plotting one feature. A \code{patchwork}
#'   object if plotting multiple features.
#' @importFrom patchwork wrap_plots
#' @importFrom stats setNames
#' @importFrom SpatialExperiment imgData getImg
#' @importMethodsFrom Matrix t
#' @export
#' @concept Spatial plotting
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
#' # Use a bounding box to zoom in
#' bbox <- c(xmin = 5500, ymin = 13500, xmax = 6000, ymax = 14000)
#' plotSpatialFeature(sfe, "nCounts", colGeometryName = "spotPoly",
#'                   annotGeometry = "myofiber_simplified",
#'                   bbox = bbox, annot_fixed = list(linewidth = 0.3))
plotSpatialFeature <- function(sfe, features, colGeometryName = 1L,
                               sample_id = "all", ncol = NULL,
                               ncol_sample = NULL,
                               annotGeometryName = NULL,
                               rowGeometryName = NULL,
                               rowGeometryFeatures = NULL,
                               annot_aes = list(), annot_fixed = list(),
                               tx_fixed = list(),
                               exprs_values = "logcounts", bbox = NULL,
                               tx_file = NULL,
                               image_id = NULL, channel = NULL, maxcell = 5e+5,
                               aes_use = c("fill", "color", "shape", "linetype"),
                               divergent = FALSE, diverge_center = NA,
                               annot_divergent = FALSE,
                               annot_diverge_center = NA,
                               size = 0.5, shape = 16, linewidth = 0,
                               linetype = 1, alpha = 1,
                               color = "black", fill = "gray80",
                               swap_rownames = NULL,
                               scattermore = FALSE, pointsize = 0,
                               bins = NULL, summary_fun = sum, hex = FALSE,
                               show_axes = FALSE, dark = FALSE,
                               palette = colorRampPalette(c("black", "white"))(255),
                               normalize_channels = FALSE,
                               ...) {
    aes_use <- match.arg(aes_use)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    values <- .get_feature_values(sfe, features, sample_id,
        colGeometryName = colGeometryName,
        exprs_values = exprs_values,
        swap_rownames = swap_rownames,
        show_symbol = !is.null(swap_rownames)
    )

    inds <- !names(values) %in% features
    if (any(inds))
        features[inds] <- rowData(sfe)[features[inds], swap_rownames]
    values <- values[,features, drop = FALSE]
    c(rowGeometryName, rowGeometryFeatures) %<-%
        .getRowGeometryFeatures(sfe, features, rowGeometryName, rowGeometryFeatures)
    .plotSpatialFeature(
        sfe, values, colGeometryName, sample_id, ncol,
        ncol_sample, annotGeometryName, annot_aes,
        annot_fixed, tx_fixed, bbox, image_id, channel, aes_use, divergent,
        diverge_center, annot_divergent,
        annot_diverge_center, size, shape, linewidth, linetype,
        alpha, color, fill, scattermore, pointsize, bins, summary_fun, hex,
        maxcell, show_axes, dark, palette, normalize_channels,
        rowGeometryName, rowGeometryFeatures, tx_file, ...
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

#' @importFrom ggplot2 geom_line theme_void theme_bw
#' @importFrom SpatialFeatureExperiment rowGeometry df2sf
.plot_graph <- function(sfe, MARGIN, sample_id, graph_name, geometry_name,
                        segment_size = 0.5, geometry_size = 0.5, ncol = NULL,
                        weights = FALSE, bbox = NULL) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (!is.null(geometry_name)) {
        gf <- switch(MARGIN,
            rowGeometry,
            colGeometry,
            annotGeometry
        )
        geometry <- gf(sfe, type = geometry_name, sample_id = sample_id)
        if (MARGIN == 2L) geometry$sample_id <- sfe$sample_id[sfe$sample_id %in% sample_id]
    } else {
        geometry <- NULL
    }
    df <- .get_graph_df(sfe, MARGIN, sample_id, graph_name, geometry)
    df <- .crop(df, bbox)
    if (!is.null(geometry)) geometry <- .crop(geometry[,"sample_id"], bbox)
    p <- ggplot()
    if (!is.null(geometry_name)) {
        p <- p +
            geom_sf(
                data = geometry, size = geometry_size, fill = NA,
                color = "gray70"
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
#' @importFrom sf st_coordinates st_centroid st_geometry st_sfc
#' @return A ggplot2 object.
#' @export
#' @concept Spatial plotting
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
plotColGraph <- function(sfe, colGraphName = 1L, colGeometryName = 1L,
                         sample_id = "all", weights = FALSE, segment_size = 0.5,
                         geometry_size = 0.5, ncol = NULL, bbox = NULL) {
    .plot_graph(sfe,
        MARGIN = 2L, sample_id = sample_id, graph_name = colGraphName,
        geometry_name = colGeometryName,
        segment_size = segment_size, geometry_size = geometry_size,
        ncol = ncol, weights = weights, bbox = bbox
    )
}

#' @rdname plotColGraph
#' @export
plotAnnotGraph <- function(sfe, annotGraphName = 1L, annotGeometryName = 1L,
                           sample_id = "all", weights = FALSE,
                           segment_size = 0.5, geometry_size = 0.5,
                           ncol = NULL, bbox = NULL) {
    .plot_graph(sfe,
        MARGIN = 3L, sample_id = sample_id, graph_name = annotGraphName,
        geometry_name = annotGeometryName,
        segment_size = segment_size, geometry_size = geometry_size,
        ncol = ncol, weights = weights, bbox = bbox
    )
}

#' Plot cell density as 2D histogram
#'
#' This function plots cell density in histological space as 2D histograms,
#' especially helpful for larger smFISH-based datasets.
#'
#' @inheritParams plotSpatialFeature
#' @param bins Number of bins. Can be a vector of length 2 to specify for x and
#'   y axes separately.
#' @param binwidth Width of bins, passed to \code{\link{geom_bin2d}} or
#'   \code{\link{geom_hex}}.
#' @param hex Logical, whether to use hexagonal bins.
#' @return A ggplot object.
#' @export
#' @concept Spatial plotting
#' @examples
#' library(SFEData)
#' sfe <- HeNSCLCData()
#' plotCellBin2D(sfe)
plotCellBin2D <- function(sfe, sample_id = "all", bins = 200, binwidth = NULL,
                          hex = FALSE, ncol = NULL, bbox = NULL) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    bin_fun <- if (hex) geom_hex else geom_bin2d
    df <- as.data.frame(spatialCoords(sfe))
    names(df) <- c("x", "y")
    df$sample_id <- colData(sfe)$sample_id
    df <- .crop(df, bbox)
    x <- y <- NULL
    p <- ggplot(df, aes(x, y)) +
        bin_fun(bins = bins, binwidth = binwidth) +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        coord_equal() +
        scale_x_continuous(expand = expansion()) +
        scale_y_continuous(expand = expansion()) +
        labs(x = NULL, y = NULL)
    if (length(sample_id) > 1L) {
        p <- p + facet_wrap(~ sample_id)
    }
    p
}

.add_image <- function(p, sfe, sample_id, image_id, channel, bbox, maxcell,
                       normalize_channels, palette) {
    if (!is.null(image_id))
        img_df <- .get_img_df(sfe, sample_id, image_id, channel, bbox, maxcell,
                              normalize_channels)
    else img_df <- NULL
    if (!is.null(image_id) && nrow(img_df)) {
        data <- NULL
        p <- p + geom_spi_rgb(data = img_df, aes(spi = data), palette = palette)
    }
    p
}

#' Plot geometries without coloring
#'
#' Different samples are plotted in separate facets. When multiple geometries
#' are plotted at the same time, they will be differentiated by color, by
#' default using the \code{dittoSeq} palette, but this can be overridden with
#' \code{scale_color_*} functions. Transcript spots of different genes are
#' differentiated by point shape if plotted, so the number of genes plotted
#' shouldn't exceed about 6 or a warning will be issued.
#'
#' @inheritParams plotSpatialFeature
#' @inheritParams SpatialFeatureExperiment::findSpatialNeighbors
#' @inheritParams plotTxBin2D
#' @param fill Logical, whether to fill polygons.
#' @param tx_alpha Transparency for transcript spots, helpful when the
#'   transcript spots are overplotting.
#' @param tx_size Point size for transcript spots.
#' @return A \code{ggplot} object.
#' @export
#' @concept Spatial plotting
#' @examples
#' library(SFEData)
#' sfe1 <- McKellarMuscleData("small")
#' sfe2 <- McKellarMuscleData("small2")
#' sfe <- cbind(sfe1, sfe2)
#' sfe <- removeEmptySpace(sfe)
#' plotGeometry(sfe, colGeometryName = "spotPoly")
#' plotGeometry(sfe, annotGeometryName = "myofiber_simplified")
plotGeometry <- function(sfe,
                         type = deprecated(), MARGIN = deprecated(),
                         colGeometryName = NULL, annotGeometryName = NULL,
                         rowGeometryName = NULL, gene = "all",
                         sample_id = "all",
                         fill = TRUE, ncol = NULL, bbox = NULL, tx_alpha = 1,
                         tx_size = 1, tx_file = NULL,
                         image_id = NULL, channel = NULL,
                         maxcell = 5e+5, show_axes = FALSE, dark = FALSE,
                         palette = colorRampPalette(c("black", "white"))(255),
                         normalize_channels = FALSE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (is_present(type))
        deprecate_warn("1.8.0", "plotGeometry(type)", details =
                           "Please use colGeometryName, annotGeometryName, or rowGeometryName instead.")
    if (is_present(MARGIN))
        deprecate_warn("1.8.0", "plotGeometry(MARGIN)", details =
                           "Please use colGeometryName, annotGeometryName, or rowGeometryName instead.")
    if (is_present(type) && is_present(MARGIN)) {
        # The old behavior
        if (MARGIN == 2L) {
            colGeometryName <- type
            annotGeometryName <- NULL
            rowGeometryName <- NULL
        }
        else if (MARGIN == 3L) {
            annotGeometryName <- type
            colGeometryName <- NULL
            rowGeometryName <- NULL
        }
    }
    if (is.null(colGeometryName) && is.null(annotGeometryName) && is.null(rowGeometryName)) {
        stop("At lease one of colGeometryName, annotGeometryName, and rowGeometryName must be specified.")
    }
    cgs <- lapply(colGeometryName, function(n) {
        out <- colGeometry(sfe, type = n, sample_id = sample_id)
        out$sample_id <- colData(sfe)[rownames(out), "sample_id"]
        out$type <- n
        out[,c("sample_id", "type", "geometry")]
    })
    ags <- lapply(annotGeometryName, function(n) {
        out <- annotGeometry(sfe, type = n, sample_id = sample_id)
        out$type <- n
        out[,c("sample_id", "type", "geometry")]
    })
    if (!is.null(rowGeometryName)) {
        if (length(rowGeometryName) > 1L) {
            rowGeometryName <- rowGeometryName[1]
            message("Only the first item in rowGeometryName is used.")
        }
        rgs <- .get_tx_df(sfe, data_dir = NULL, tech = NULL, file = NULL,
                          sample_id = sample_id, spatialCoordsNames = c("X", "Y"),
                          gene_col = "gene", bbox = bbox, gene = gene,
                          return_sf = TRUE, rowGeometryName = rowGeometryName,
                          geoparquet_file = tx_file)
    } else rgs <- NULL
    df <- do.call(rbind, c(cgs, ags))
    p <- ggplot()
    p <- .add_image(p, sfe, sample_id, image_id, channel, bbox, maxcell,
                    normalize_channels, palette)
    if (dark) color <- "gray90" else color <- "black"
    if (!is.null(rowGeometryName)) {
        if (!identical(gene, "all") && length(unique(rgs$gene)) > 1L)
            p <- p + geom_sf(data = rgs, aes(shape = gene), color = color,
                             alpha = tx_alpha, size = tx_size,
                             show.legend = "point")
        else p <- p + geom_sf(data = rgs, shape = 4, color = color,
                              alpha = tx_alpha, size = tx_size,
                              show.legend = "point")
    }
    if (!is.null(df)) {
        df <- .crop(df, bbox)
        if (fill) {
            if (dark) fill_col <- "darkblue" else fill_col <- "gray90"
        } else fill_col <- NA
        if (length(unique(df$type)) > 1L) {
            p <- p + geom_sf(data = df, aes(color = type), fill = NA) +
                scale_color_manual(values = ditto_colors)
        } else {
            p <- p + geom_sf(data = df, fill = fill_col, color = color)
        }
    }
    if (length(sample_id) > 1L) {
        p <- p + facet_wrap(~ sample_id, ncol = ncol)
    }
    if (dark) p <- p + .dark_theme(show_axes)
    else if (show_axes) p <- p + theme_bw() else p <- p + theme_void()
    p
}

#' Show image without plotting geometries
#'
#' This function plots the images in SFE objects without plotting geometries.
#' When showing axes, the numbers are coordinates within the image itself and
#' have the same units as the spatial extent, but are not the actual spatial
#' extent when plotting multiple samples to avoid excessive empty space.
#'
#' @inheritParams plotGeometry
#' @param image_id ID of the image(s) to plot. If \code{NULL}, then the first
#'   image present is plotted. Can be a vector of IDs to use different grayscale
#'   images for different channels. The vector can be named ('r', 'g', 'b'), to
#'   assign channels to images. The vector must be named if it's length 2.
#' @return A \code{ggplot} object.
#' @concept Spatial plotting
#' @export
#' @examples
#' library(SFEData)
#' library(SpatialFeatureExperiment)
#' fn <- XeniumOutput("v2", file_path = "xenium_example")
#' # Weird RBioFormats null pointer error the first time it's run
#' try(sfe <- readXenium(fn))
#' sfe <- readXenium(fn)
#' # Plot one channel
#' plotImage(sfe, image_id = "morphology_focus", channel = 1L)
#' plotImage(sfe, image_id = "morphology_focus", channel = 1L, show_axes = TRUE, dark = TRUE)
#' # Colorize based on different channels
#' plotImage(sfe, image_id = "morphology_focus", channel = c(2,4,1), show_axes = TRUE, dark = TRUE)
#' unlink("xenium_example", recursive = TRUE)
plotImage <- function(sfe, sample_id = "all", image_id = NULL, channel = NULL,
                      ncol = NULL, bbox = NULL,
                      maxcell = 5e+5, show_axes = FALSE, dark = FALSE,
                      palette = colorRampPalette(c("black", "white"))(255),
                      normalize_channels = FALSE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    img_df <- .get_img_df(sfe, sample_id, image_id, channel, bbox, maxcell,
                          normalize_channels)
    data <- NULL
    img_df$bbox <- lapply(img_df$data, function(x) {
        out <- ext(x)
        mins <- out[c("xmin", "ymin")]
        out <- st_bbox(out) |> st_as_sfc()
        if (length(sample_id) > 1L) {
            out <- out - mins
        }
        out
    })
    img_df$bbox <- st_sfc(unlist(img_df$bbox, recursive = FALSE))
    p <- ggplot() +
        geom_sf(data = img_df, aes(geometry = bbox), fill = NA, linewidth = 0) +
        geom_spi_rgb(data = img_df, aes(spi = data), palette = palette)
    if (length(sample_id) > 1L) {
        p <- p + facet_wrap(~ sample_id, ncol = ncol)
    }
    if (dark) p <- p + .dark_theme(show_axes)
    else if (show_axes) p <- p + theme_bw() else p <- p + theme_void()
    p
}
