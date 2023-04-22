# From tidyterra and ggspavis::plotVisium
resample_spat <- function(r, maxcell = 50000) {
    if (terra::ncell(r) > 1.1 * maxcell) {
        r <- terra::spatSample(r, maxcell,
                               as.raster = TRUE,
                               method = "regular"
        )
    }
    return(r)
}

# Instead of one SpatRaster, data should be a data frame with a list column for
# the images as SpatRasterImage and another column sample_id
geom_spi_rgb <- function(mapping = NULL,
                         data = NULL,
                         interpolate = TRUE,
                         maxcell = 500000,
                         max_col_value = 255,
                         ...) {

    layer_spatrast <- ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = StatSpiRGB,
        geom = GeomSpiRGB,
        position = "identity",
        inherit.aes = FALSE,
        show.legend = FALSE,
        params = list(
            na.rm = TRUE,
            # Extra params
            maxcell = maxcell,
            interpolate = interpolate,
            max_col_value = max_col_value,
            ...
        )
    )

    layer_spatrast
}

# Stats----
StatSpiRGB <- ggplot2::ggproto(
    "StatSpiRGB",
    ggplot2::Stat,
    required_aes = "spi",
    extra_params = c("maxcell", "max_col_value", "na.rm"),
    compute_layer = function(self, data, params, layout) {
        ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_layer(
            data,
            params, layout
        )
    },
    compute_group = function(data, scales, coord, params,
                             maxcell = 5e+5) {
        # Can only plot one image per facet
        data <- data[1,]
        # Extract raster from group
        rast <- imgRaster(data$spi[[1]])
        if (length(names(rast)) == 3L) {
            names(rast) <- c("r", "g", "b")
            # Remove RGB settings, better plot without it
            terra::RGB(rast) <- NULL
        }
        rast <- resample_spat(rast, maxcell)
        data$spi <- I(list(rast))
        data
    }
)


# Geom----
#' @importFrom terra ext
#' @importFrom grDevices as.raster
#' @importFrom grid rasterGrob
# Based on geom_raster ggplot2
GeomSpiRGB <- ggplot2::ggproto(
    "GeomSpiRGB",
    ggplot2::GeomCustomAnn,
    required_aes = c("spi"),
    draw_panel = function(self, data, panel_params, coord, interpolate = FALSE) {
        bbox <- as.vector(ext(data$spi[[1]]))
        if (length(names(data$spi[[1]])) == 3L)
            raster <- terra::as.array(data$spi[[1]]) |> as.raster(max = 255)
        else
            raster <- terra::as.array(data$spi[[1]])[,,1] |> as.raster(max = 255)
        ggplot2::ggproto_parent(ggplot2::GeomCustomAnn, self)$draw_panel(
            panel_params = panel_params, coord = coord,
            xmin = bbox["xmin"], xmax = bbox["xmax"],
            ymin = bbox["ymin"], ymax = bbox["ymax"],
            grob = rasterGrob(raster, interpolate = interpolate)
        )
    }
)

geom_spi <- function(mapping = NULL,
                     data = NULL,
                     na.rm = TRUE,
                     show.legend = NA,
                     inherit.aes = FALSE,
                     interpolate = FALSE,
                     maxcell = 500000,
                     ...) {

    # Create layer
    ggplot2::layer(
        data = data,
        mapping = modifyList(mapping, aes(fill = ggplot2::after_stat(.data$value))),
        stat = StatSpi,
        geom = ggplot2::GeomRaster,
        position = "identity",
        inherit.aes = inherit.aes,
        show.legend = show.legend,
        params = list(
            na.rm = na.rm,
            # Extra params
            maxcell = maxcell,
            interpolate = interpolate,
            ...
        )
    )
}


StatSpi <- ggplot2::ggproto(
    "StatSpi",
    ggplot2::Stat,
    required_aes = "spi",
    extra_params = c("maxcell", "na.rm"),
    compute_layer = function(self, data, params, layout) {
        ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_layer(
            data,
            params, layout
        )
    },
    compute_group = function(data, scales, coord, params, maxcell = 5e+5) {
        # Can only plot one image per facet
        data <- data[1,]
        # Extract raster from group
        rast <- imgRaster(data$spi[[1]])
        if (length(names(rast)) > 1L) {
            warning("Only the first layer of the image is plotted.")
            rast <- rast[[1]]
        }
        rast <- resample_spat(rast, maxcell)
        out <- terra::as.data.frame(rast, xy = TRUE)
        names(out)[3] <- "value"
        out$sample_id <- data$sample_id
        out
    }
)
