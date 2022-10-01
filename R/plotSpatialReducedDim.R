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
#' @return Same as in \code{\link{plotSpatialFeature}}. A \code{ggplot2} object
#'   if plotting one component. A \code{patchwork} object if plotting multiple
#'   components.
#' @export
#' @seealso scater::plotReducedDim
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- logNormCounts(sfe)
#' sfe <- runPCA(sfe, ncomponents = 2)
#' spatialReducedDim(sfe, "PCA", 2, "spotPoly")
spatialReducedDim <- function(sfe, dimred, ncomponents, colGeometryName = 1L,
                              sample_id = NULL, ncol = NULL, ncol_sample = NULL,
                              aes_use = c("fill", "color", "shape", "linetype"),
                              divergent = FALSE, diverge_center = NULL,
                              size = 0, shape = 16, linetype = 1, alpha = 1,
                              color = NA, fill = "gray80", show_symbol = TRUE,
                              ...) {
    aes_use <- match.arg(aes_use)
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    if (length(ncomponents) == 1L) {
        dims_use <- seq_len(ncomponents)
    } else {
        dims_use <- ncomponents
    }
    sample_ind <- colData(sfe)$sample_id %in% sample_id
    values <- as.data.frame(reducedDim(sfe, dimred)[sample_ind,dims_use])
    .plotSpatialFeature(sfe, values, colGeometryName, sample_id, ncol,
                        ncol_sample, annotGeometryName = NULL,
                        annot_aes = list(), annot_fixed = list(), aes_use,
                        divergent, diverge_center, annot_divergent = FALSE,
                        annot_diverge_center = NULL, size, shape, linetype,
                        alpha, color, fill, show_symbol,...)
}
