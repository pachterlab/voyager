#' Convert listw into sparse adjacency matrix
#'
#' Edge weights are used in the adjacency matrix. Because most elements of the
#' matrix are 0, using sparse matrix greatly reduces memory use.
#'
#' @param listw A \code{listw} object for spatial neighborhood graph.
#' @return A sparse \code{dgCMatrix}, whose row represents each cell or spot and
#'   whose columns represent the neighbors. The matrix does not have to be
#'   symmetric.
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
    sparseMatrix(i = i, j = j, x = x)
}
