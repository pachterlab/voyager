# 6. Find high-high, low-low, high-low, low-high regions from moran.plot and localmoran

#' Moran scatterplot
#'
#' To make it easier to compute the Moran scatterplot for a bunch of genes at
#' once. And to reduce the boilerplate for getting the gene expression values.
#' Burning question: how shall I compare these across genes? Also to organize
#' the results, which was one of the motivations of SFE.
#'
#' @inheritParams calculateMoransI
#' @param ... Other arguments passed to \code{\link{moran.plot}}.
#' @name calculateMoranPlot
#' @importFrom spdep moran.plot
#' @return For \code{calculateMoranPlot}, a list of data frames that are the
#'   output of \code{\link{moran.plot}}. The names of the list are the names of
#'   the features. For \code{runMoranPlot}, the list is added to a column in
#'   \code{rowData(x)}. For the colData, colGeometry, and annotGeometry
#'   versions, the results are added to an attribute of the data frame of
#'   interest called \code{featureData}, in a manner analogous to
#'   \code{rowData}. The plot is not made.
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute Moran plot for vector or matrix
#' calculateMoranPlot(colData(sfe)$nCounts, listw = colGraph(sfe, "visium"))
#' # Add results to rowData, features are genes
#' sfe <- runMoranPlot(sfe, features = rownames(sfe)[1], exprs_values = "counts")
#' rowData(sfe)
#' # Specifically for colData
#' sfe <- colDataMoranPlot(sfe, "nCounts")
#' attr(colData(sfe), "featureData")
NULL

#' @rdname calculateMoranPlot
#' @export
setMethod("calculateMoranPlot", "ANY",
          function(x, listw, BPPARAM = SerialParam(),
                   zero.policy = NULL, ...) {
            .calc_univar_autocorr(x, listw, fun = moran.plot, BPPARAM, returnDF = FALSE,
                                  plot = FALSE, return_df = TRUE,
                                  zero.policy = zero.policy, ...)
          })

#' @rdname calculateMoranPlot
#' @export
setMethod("calculateMoranPlot", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateMoranPlot, returnDF = FALSE))

.MoranPlot2df <- function(out, name) {
  out_df <- DataFrame(res = I(out))
  names(out_df) <- name
  out_df
}

#' @rdname calculateMoranPlot
#' @export
colDataMoranPlot <- .coldata_univar_fun(calculateMoranPlot, .MoranPlot2df, "MoranPlot")

#' @rdname calculateMoranPlot
#' @export
colGeometryMoranPlot <- .colgeom_univar_fun(calculateMoranPlot, .MoranPlot2df, "MoranPlot")

#' @rdname calculateMoranPlot
#' @export
annotGeometryMoranPlot <- .annotgeom_univar_fun(calculateMoranPlot, .MoranPlot2df, "MoranPlot")

#' @rdname calculateMoranPlot
#' @export
runMoranPlot <- function(sfe, features, colGraphName = 1L, sample_id = NULL,
                         exprs_values = "logcounts", BPPARAM = SerialParam(),
                         zero.policy = NULL, name = "MoranPlot", ...) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  for (s in sample_id) {
    out <- calculateMoranPlot(sfe, features, colGraphName, s, exprs_values,
                              BPPARAM, zero.policy, ...)
    out <- .MoranPlot2df(out, name)
    features <- .symbol2id(sfe, features)
    rownames(out) <- features
    out <- .add_name_sample_id(out, s)
    rowData(sfe)[features, names(out)] <- out
  }
  sfe
}

#' Find clusters on the Moran plot
#'
#' The Moran plot plots the value at each location on the x axis, and the
#' average of the neighbors of each locations on the y axis. Sometimes clusters
#' can be seen on the Moran plot, indicating different types of neighborhoods.
#'
#' @inheritParams bluster::clusterRows
#' @inheritParams calculateMoranPlot
#' @inheritParams calculateMoransI
#' @inheritParams moranPlot
#' @param sfe A \code{SpatialFeatureExperiment} object with Moran plot computed
#'   for the feature of interest. If the Moran plot for that feature has not
#'   been computed for that feature in this sample_id, it will be calculated and
#'   stored in \code{rowData}. See \code{\link{calculateMoranPlot}}.
#' @param colGeometryName Name of colGeometry from which to look for features.
#' @param annotGeometryName Name of annotGeometry from which to look for
#'   features.
#' @param features Features whose Moran plot are to be cluster. Features whose
#'   Moran plots have not been computed will be skipped, with a warning.
#' @return A \code{DataFrame} each column of which is a factor for cluster
#'   membership of each feature. The column names are the features.
#' @importFrom bluster clusterRows
#' @importFrom methods as
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' library(bluster)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute Moran plot for vector or matrix
#' calculateMoranPlot(colData(sfe)$nCounts, listw = colGraph(sfe, "visium"))
#' # Add results to rowData, features are genes
#' sfe <- runMoranPlot(sfe, features = rownames(sfe)[1], exprs_values = "counts")
#' sfe <- colDataMoranPlot(sfe, "nCounts")
#' clusts <- clusterMoranPlot(sfe, c(rownames(sfe)[1], "nCounts"),
#'                            BLUSPARAM = KmeansParam(2))
clusterMoranPlot <- function(sfe, features, BLUSPARAM, sample_id = NULL,
                             name = "MoranPlot",
                             colGeometryName = NULL,
                             annotGeometryName = NULL, show_symbol = TRUE) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  use_col <- is.null(annotGeometryName) || !is.null(colGeometryName)
  out <- lapply(sample_id, function(s) {
    colname_use <- paste(name, s, sep = "_")
    mps <- .get_feature_metadata(sfe, features, name, s, colGeometryName,
                                 annotGeometryName, show_symbol)
    o <- lapply(mps, function(mp) clusterRows(mp[,c("x", "wx")], BLUSPARAM))
    o <- as(o, "DataFrame")
    o$sample_id <- s
    if (use_col)
      row.names(o) <- colnames(sfe)[colData(sfe)$sample_id == s]
    # What if some features don't have the Moran Plot computed
    features_absent <- setdiff(features, names(o))
    if (length(features_absent) && length(sample_id) > 1L) {
      for (f in features_absent) {
        o[[f]] <- NA
      }
      o <- o[,c("sample_id", features)]
    }
    o
  })
  if (length(sample_id) > 1L) {
    out <- do.call(rbind, out)
    for (f in features) {
      if (all(is.na(out[[f]]))) out[[f]] <- NULL
      else {
        out[[f]] <- factor(out[[f]], levels = seq_len(max(as.integer(out[[f]]))))
      }
    }
  } else out <- out[[1]]
  out
}
