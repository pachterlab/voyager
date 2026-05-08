#' Compute the bounds of Moran's I given spatial neighborhood graph
#'
#' Values Moran's I can take depends on the spatial neighborhood graph. The
#' bounds of Moran's I given the graph, C, are given by the minimum and maximum
#' eigenvalues of the double centered -- i.e. subtracting column means and row
#' means -- adjacency matrix \eqn{(I - \mathbb{11}^T/n)C(I - \mathbb{11}^T/n)}{(I
#' - 11^T/n)C(I - 11^T/n)}, where \eqn{\mathbb 1}{1} is a vector of all 1's.
#' This implementation follows the implementation in \code{adespatial} and uses
#' the \code{RSpectra} package to more quickly find only the minimum and maximum
#' eigenvalues without performing unnecessary work to find the full spectrum as
#' done in base R's \code{\link{eigen}}.
#'
#' @param listw A listw object for the spatial neighborhood graph.
#' @return A numeric vector of minimum and maximum Moran's I given the spatial
#'   neighborhood graph.
#' @note
#' After double centering, the adjacency matrix is no longer sparse, so this
#' function can take up a lot of memory for larger datasets.
#' @export
#' @importFrom RSpectra eigs_sym
#' @importFrom spdep listw2mat
#' @importFrom SpatialFeatureExperiment multi_listw2sparse
#' @concept Spatial statistics
#' @references de Jong, P., Sprenger, C., & van Veen, F. (1984). On extreme values of Moran's I and Geary's C. Geographical Analysis, 16(1), 17-24.
#' @examples
#' # example code
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' g <- findVisiumGraph(sfe)
#' moranBounds(g)
moranBounds <- function(listw) {
    W <- listw2mat(listw)
    W <- (W + t(W)) / 2
    rs <- nrow(W)/sum(W)
    # Double center
    row_means <- rowMeans(W)
    col_means <- colMeans(W) - mean(row_means)
    W <- sweep(W, 2, row_means)
    W <- sweep(W, 1, col_means)
    out <- rs * sort(eigs_sym(W, k = 2, which = "BE", opts = list(retvec = FALSE))$values)
    names(out) <- c("Imin", "Imax")
    out
}
