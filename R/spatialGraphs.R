# Compute and manage spatial neighborhood graphs
# Giotto can do triangulation, knn, but spdep has more methods
# I don't think Seurat uses the spatial graph
# 1. Wrap the neighborhood methods in spdep for SFE objects
# 2. Remember polygon contiguity. Also try st_is_within_distance
# 3. For Visium, see if vi2mrc works, and if it does, whether it's faster than using distance
# 4. Check distribution of distance between spots in units associated with the SFE object.
# 5. Diagnostic functions for graphs such as cardinality and number/proportion of singletons
# To be made easy to call with the SFE object.
# These internal wrapper functions return nb, which will later be converted to listw
.tri2nb_sfe <- function(x) tri2nb(spatialCoords(x))
.knn_sfe <- function(x, k = 1, use_kd_tree = TRUE)
  knn2nb(knearneigh(spatialCoords(x), k = k, use_kd_tree = use_kd_tree,
                    longlat = FALSE))
.dnn_sfe <- function(x, d1, d2, use_kd_tree = TRUE)
  dnearneigh(spatialCoords(x), d1, d2, longlat = FALSE,
             use_kd_tree = use_kd_tree)
.g2nb_sfe <- function(x, fun, nnmult = 3, sym = FALSE) {
  # Either gabrielneigh or relativeneigh
  g <- fun(spatialCoords(x), nnmult)
  graph2nb(g, sym = sym)
}
.gabriel_sfe <- function(x, nnmult = 3, sym = FALSE) .g2nb_sfe(x, gabrielneigh, nnmult, sym)
.relative_sfe <- function(x, nnmult = 3, sym = FALSE) .g2nb_sfe(x, relativeneigh, nnmult, sym)
.soi_sfe <- function(x, quadsegs = 10, sym = FALSE) {
  coords <- spatialCoords(x)
  g <- soi.graph(tri2nb(coords), coords, quadsegs)
  graph2nb(g, sym = sym)
}

#' Find spatial neighborhood graph
#'
#' This function wraps all spatial neighborhood graphs implemented in the
#' package \code{spdep} for the \code{SpatialFeatureExperiment} (SFE) class, to
#' find spatial neighborhood graphs for the entities represented by columns of
#' the gene count matrix in the SFE object. Results are stored as \code{listw}
#' objects in the \code{spatialGraphs} field of the SFE object, as \code{listw}
#' is used in many methods that facilitate the spatial neighborhood graph in the
#' \code{spdep}, \code{spatialreg}, and \code{adespatial}. The edge weights of
#' the graph in the \code{listw} object are by default style W (see
#' \code{\link{nb2listw}}) and the unweighted neighbor list is in the
#' \code{neighbours} field of the \code{listw} object.
#'
#' @inheritParams spdep::nb2listw
#' @param x A \code{\link{SpatialFeatureExperiment}} object.
#' @param method Name of function in the package \code{spdep} to use to find the
#'   spatial neighborhood graph.
#' @param name Name of the graph to store in the \code{spatialGraphs} of the
#'   output. Defaults to the name of the function in \code{method}.
#' @param ... Extra arguments passed to the \code{spdep} function stated in the
#'   \code{method} argument, such as \code{k}, \code{use_kd_tree}, \code{d1},
#'   \code{d2}, \code{nnmult}, \code{sym}, and \code{quadsegs}. Note that any
#'   arguments about using longitude and latitude, which are irrelevant, are
#'   ignored. The \code{longlat} argument is hard coded to \code{FALSE}.
#' @return A \code{SpatialFeatureExperiment} object with spatial neighborhood
#'   graph.
#' @importFrom spdep tri2nb knearneigh dnearneigh gabrielneigh relativeneigh
#' soi.graph knn2nb graph2nb nb2listw
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SpatialFeatureExperiment spatialGraph
#' @export
setMethod("findSpatialNeighbors", "SpatialFeatureExperiment",
          function(x, method = c("tri2nb", "knearneigh", "dnearneigh",
                                 "gabrielneigh", "relativeneigh", "soi.graph"),
                   name = method, glist = NULL, style = "W", zero.policy = NULL,
                   ...) {
            method <- match.arg(method)
            extra_args_use <- switch (method,
              tri2nb = NULL,
              knearneigh = c("k", "use_kd_tree"),
              dnearneigh = c("d1", "d2", "use_kd_tree"),
              gabrielneigh = c("nnmult", "sym"),
              relativeneigh = c("nnmult", "sym"),
              soi.graph = c("quadsegs", "sym")
            )
            args <- list(...)
            args <- args[names(args) %in% extra_args_use]
            fun_use <- switch (method,
              tri2nb = .tri2nb_sfe,
              knearneigh = .knn_sfe,
              dnearneigh = .dnn_sfe,
              gabrielneigh = .gabriel_sfe,
              relativeneigh = .relative_sfe,
              soi.graph = .soi_sfe
            )
            nb_out <- do.call(fun_use, c(x = x, args))
            listw_use <- nb2listw(nb_out, glist, style, zero.policy)
            spatialGraph(x, name) <- listw_use
            x
          })

#' Find spatial neighborhood graphs for Visium spots
#'
#' Visium spots are arranged in a hexagonal grid. This function uses the known
#' locations of the Visium barcodes to construct a neighborhood graph, so
#' adjacent spots are connected by edges. Since the known rows and columns of
#' the spots are used, the unit the spot centroid coordinates are in does not
#' matter.
#'
#' @inheritParams spdep::nb2listw
#' @param x A \code{SpatialFeatureExperiment} object with Visium data. Column
#'   names of the gene count matrix must be Visium barcodes, which may have a
#'   numeric suffix to distinguish between samples (e.g. "AAACAACGAATAGTTC-1").
#' @param sample_id String specifying the sample, as in \code{sample_id} in
#'   \code{colData}, for which the graph is to be constructed.
#' @return A \code{listw} object for the neighborhood graph. The node index will
#' be the order of barcodes within the sample in the gene count matrix.
#' @importFrom spdep dnearneigh nb2listw
#' @export
findVisiumGraph <- function(x, sample_id, style = "W", zero.policy = NULL) {
  bcs_use <- colnames(x)[colData(x)$sample_id == sample_id]
  bcs_use2 <- sub("-\\d+$", "", bcs_use)
  coords_use <- visium_row_col[match(bcs_use2, visium_row_col), c("col", "row")]
  # So adjacent spots are equidistant
  coords_use$row <- coords_use$row * sqrt(3)
  g <- dnearneigh(as.matrix(coords_use), d1 = 1.9, d2 = 2.1, row.names = bcs_use)
  out <- nb2listw(g, style = style, zero.policy = zero.policy)
  return(out)
}
