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
                         palette = colorRampPalette(c("black", "white"))(255),
                         ...) {

    layer_spatrast <- ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = ggplot2::StatIdentity,
        geom = GeomSpiRGB,
        position = "identity",
        inherit.aes = FALSE,
        show.legend = FALSE,
        params = list(
            na.rm = TRUE,
            # Extra params
            interpolate = interpolate,
            palette = palette
        )
    )

    layer_spatrast
}

# Geom----
#' @importFrom terra ext as.raster ncell nlyr
#' @importFrom grDevices colorRampPalette
#' @importFrom grid rasterGrob
# Based on geom_raster ggplot2
GeomSpiRGB <- ggplot2::ggproto(
    "GeomSpiRGB",
    ggplot2::GeomCustomAnn,
    required_aes = c("spi"),
    default_aes = aes(channel = NA),
    draw_panel = function(self, data, panel_params, coord,
                          palette=colorRampPalette(c("black", "white"))(255),
                          interpolate = FALSE) {
        img <- data$spi[[1]]
        bbox <- ext(img)
        max_value <- terra::minmax(img, compute = TRUE)[2,] |> max()
        # Would be nice to allow other color maps
        if (nlyr(img) == 3L)
            raster <- terra::as.array(img) |> as.raster(max = max_value)
        else
            raster <- terra::as.raster(img, maxcell = Inf, col = palette)
        ggplot2::ggproto_parent(ggplot2::GeomCustomAnn, self)$draw_panel(
            panel_params = panel_params, coord = coord,
            xmin = bbox["xmin"], xmax = bbox["xmax"],
            ymin = bbox["ymin"], ymax = bbox["ymax"],
            grob = rasterGrob(raster, interpolate = interpolate)
        )
    }
)
