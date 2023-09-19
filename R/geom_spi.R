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
        img <- data$spi[[1]]
        bbox <- as.vector(ext(img))
        # Check if 8 bit or 16 bit
        mm <- terra::minmax(img, compute = TRUE)
        if (mm[2] < 256) max_use <- 256 else max_use <- 2^16
        if (length(names(img)) == 3L)
            raster <- terra::as.array(img) |> as.raster(max = max_use)
        else
            raster <- terra::as.array(img)[,,1] |> as.raster(max = max_use)
        ggplot2::ggproto_parent(ggplot2::GeomCustomAnn, self)$draw_panel(
            panel_params = panel_params, coord = coord,
            xmin = bbox["xmin"], xmax = bbox["xmax"],
            ymin = bbox["ymin"], ymax = bbox["ymax"],
            grob = rasterGrob(raster, interpolate = interpolate)
        )
    }
)
