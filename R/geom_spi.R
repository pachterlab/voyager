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
            ...
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
                          col=colorRampPalette(c("black", "white"))(255),
                          interpolate = FALSE) {
        img <- data$spi[[1]]
        bbox <- as.vector(ext(img))
        # Check if 8 bit or 16 bit
        mm <- terra::minmax(img, compute = TRUE)
        #if (mm[2] < 256) max_use <- 256 else max_use <- 2^16
        # Would be nice to allow other color maps
        if (length(nlyr(img)) == 3L)
            raster <- terra::as.array(img) |> as.raster(max = mm[2])
        else
            raster <- terra::as.raster(img, maxcell = Inf, col = col)
        ggplot2::ggproto_parent(ggplot2::GeomCustomAnn, self)$draw_panel(
            panel_params = panel_params, coord = coord,
            xmin = bbox["xmin"], xmax = bbox["xmax"],
            ymin = bbox["ymin"], ymax = bbox["ymax"],
            grob = rasterGrob(raster, interpolate = interpolate)
        )
    }
)
