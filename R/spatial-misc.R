#' Convert listw into sparse adjacency matrix
#'
#' Edge weights are used in the adjacency matrix. Because most elements of the
#' matrix are 0, using sparse matrix greatly reduces memory use.
#'
#' @param listw A \code{listw} object for spatial neighborhood graph.
#' @return A sparse \code{dgCMatrix}, whose row represents each cell or spot and
#'   whose columns represent the neighbors. The matrix does not have to be
#'   symmetric. If \code{region.id} is present in the \code{listw} object, then
#'   it will be the row and column names of the output matrix.
#' @export
#' @importFrom Matrix sparseMatrix
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' g <- findVisiumGraph(sfe)
#' mat <- listw2sparse(g)
listw2sparse <- function(listw) {
    i <- rep(seq_along(listw$neighbours), times = card(listw$neighbours))
    j <- unlist(listw$neighbours)
    x <- unlist(listw$weights)
    n <- length(listw$neighbours)
    region_id <- attr(listw$neighbours, "region.id")
    sparseMatrix(i = i, j = j, x = x, dims = rep(n, 2),
                 dimnames = list(region_id, region_id))
}

#' Convert multiple listw graphs into a single sparse adjacency matrix
#'
#' Each sample in the SFE object has a separate spatial neighborhood graph.
#' Spatial analyses performed jointly on multiple samples require a combined
#' spatial neighborhood graph from the different samples, where the different
#' samples would be disconnected components of the graph. This combined
#' adjacency matrix can be used in MULTISPATI PCA.
#'
#' @param listws A list of \code{listw} objects.
#' @return A sparse \code{dgCMatrix} of the combined spatial neighborhood graph,
#' with the original spatial neighborhood graphs of the samples on the diagonal.
#' When the input is an SFE object, the rows and columns will match the column
#' names of the SFE object.
#' @export
#' @examples
#' # example code
#'
multi_listw2sparse <- function(listws) {
    slices <- list()
    n <- length(listws)
    mats <- lapply(listws, listw2sparse)
    ncells <- vapply(mats, nrow, FUN.VALUE = integer(1))
    region_ids <- lapply(listws, function(l) attr(l$neighbours, "region.id"))
    tot <- sum(ncells)
    prev <- 0
    next_n <- tot
    prev_inds <- 0
    next_inds <- 1
    for (i in seq_along(listws)) {
        n_curr <- ncells[i]
        next_n <- next_n - n_curr
        next_inds <- next_inds + 1
        if (prev > 0) {
            prev_m <- .empty_dgc(nrow = prev, ncol = n_curr)
            rownames(prev_m) <- unlist(region_ids[seq_len(prev_inds)])
            o <- rbind(prev_m, mats[[i]])
        } else o <- mats[[i]]
        if (next_n > 0) {
            next_m <- .empty_dgc(nrow = next_n, ncol = n_curr)
            rownames(next_m) <- unlist(region_ids[seq(next_inds, n, by = 1)])
            o <- rbind(o, next_m)
        }
        slices[[i]] <- o
        prev <- prev + n_curr
        prev_inds <- prev_inds + 1
    }
    do.call(cbind, slices)
}

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
#' @inheritParams listw2sparse
#' @return A numeric vector of minimum and maximum Moran's I given the spatial
#'   neighborhood graph.
#' @note
#' After double centering, the adjacency matrix is no longer sparse, so this
#' function can take up a lot of memory for larger datasets.
#' @export
#' @importFrom RSpectra eigs_sym
#' @importFrom spdep listw2mat
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
