#' Plot dimension reduction components in space
#'
#' Such as plotting the value of projection of gene expression of each cell to a
#' principal component in space. At present, this function does not work for the
#' 3D array of geographically weighted PCA (GWPCA), but a future version will
#' deal with GWPCA results.
#'
#' @inheritParams plotSpatialFeature
#' @param dimred A string or integer scalar indicating the reduced dimension
#'   result in \code{reducedDims(sfe)} to plot.
#' @param ncomponents A numeric scalar indicating the number of dimensions to
#'   plot, starting from the first dimension. Alternatively, a numeric vector
#'   specifying the dimensions to be plotted.
#' @param components A numeric scalar or vector specifying which dimensions to
#'   be plotted. Use this instead of \code{ncomponents} when plotting only one
#'   dimension.
#' @return Same as in \code{\link{plotSpatialFeature}}. A \code{ggplot2} object
#'   if plotting one component. A \code{patchwork} object if plotting multiple
#'   components.
#' @export
#' @seealso scater::plotReducedDim
#' @concept Spatial plotting
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- logNormCounts(sfe)
#' sfe <- runPCA(sfe, ncomponents = 2)
#' spatialReducedDim(sfe, "PCA", ncomponents = 2, "spotPoly",
#'     annotGeometryName = "tissueBoundary",
#'     divergent = TRUE, diverge_center = 0
#' )
#' # Basically PC1 separates spots not on tissue from those on tissue.
spatialReducedDim <- function(sfe, dimred, ncomponents = NULL,
                              components = ncomponents, colGeometryName = 1L,
                              sample_id = "all", ncol = NULL, ncol_sample = NULL,
                              annotGeometryName = NULL,
                              annot_aes = list(), annot_fixed = list(),
                              exprs_values = "logcounts", bbox = NULL,
                              image_id = NULL, channel = NULL, maxcell = 5e+5,
                              aes_use = c("fill", "color", "shape", "linetype"),
                              divergent = FALSE, diverge_center = NULL,
                              annot_divergent = FALSE,
                              annot_diverge_center = NULL,
                              size = 0, shape = 16, linewidth = 0,
                              linetype = 1, alpha = 1,
                              color = NA, fill = "gray80", scattermore = FALSE,
                              pointsize = 0, bins = NULL, summary_fun = sum,
                              hex = FALSE, show_axes = FALSE, dark = FALSE,
                              palette = colorRampPalette(c("black", "white"))(255),
                              normalize_channels = FALSE, ...) {
    aes_use <- match.arg(aes_use)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (length(ncomponents) == 1L) {
        dims_use <- seq_len(ncomponents)
    } else if (length(ncomponents) > 1L) {
        dims_use <- ncomponents
    } else dims_use <- components
    sample_ind <- colData(sfe)$sample_id %in% sample_id
    values <- as.data.frame(reducedDim(sfe, dimred)[sample_ind, dims_use, drop = FALSE])
    out <- .plotSpatialFeature(sfe, values, colGeometryName, sample_id, ncol,
        ncol_sample, annotGeometryName,
        annot_aes, annot_fixed, bbox, image_id, aes_use,
        divergent, diverge_center, annot_divergent,
        annot_diverge_center, size, shape, linewidth, linetype,
        alpha, color, fill,
        show_symbol = FALSE, scattermore = scattermore, pointsize = pointsize,
        bins = bins, summary_fun = summary_fun, hex = hex, maxcell = maxcell,
        channel = channel, show_axes = show_axes, dark = dark, palette = palette,
        normalize_channels = normalize_channels, ...
    )
    if (is(out, "patchwork")) {
        out <- out + plot_annotation(title = dimred)
    } else {
        out <- out + ggtitle(dimred)
    }
    out
}
