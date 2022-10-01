#' Find clusters of correlogram patterns
#'
#' Cluster the correlograms to find patterns in length scales of spatial
#' autocorrelation. All the correlograms clustered must be computed with the
#' same method and have the same number of lags.
#'
#' @inheritParams clusterMoranPlot
#' @inheritParams calculateUnivariate
#' @inheritParams plotCorrelogram
#' @param sfe A \code{SpatialFeatureExperiment} object with correlograms
#' computed for features of interest.
#' @param features Features whose correlograms to cluster.
#' @return A \code{DataFrame} with 3 columns: \code{feature} for the features,
#' \code{cluster} a factor for cluster membership of the features within each
#' sample, and \code{sample_id} for the sample.
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SFEData)
#' library(bluster)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' inds <- c(1,3,4,5)
#' sfe <- runUnivariate(sfe, type = "sp.correlogram",
#'                      features = rownames(sfe)[inds],
#'                      exprs_values = "counts", order = 5)
#' clust <- clusterCorrelograms(sfe, features = rownames(sfe)[inds],
#'                              BLUSPARAM = KmeansParam(2))
clusterCorrelograms <- function(sfe, features, BLUSPARAM, sample_id = NULL,
                                method = "I",
                                colGeometryName = NULL,
                                annotGeometryName = NULL, show_symbol = TRUE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    name <- paste("sp.correlogram", method, sep = "_")
    out <- lapply(sample_id, function(s) {
        ress <- .get_feature_metadata(sfe, features, name, s, colGeometryName,
                                      annotGeometryName, show_symbol = show_symbol)
        if (method %in% c("I", "C")) {
            # First column is the metric, second column expectation, third is variance
            ress <- lapply(ress, function(r) r[,1])
        }
        res_mat <- do.call(rbind, ress)
        rownames(res_mat) <- names(ress)
        clus <- clusterRows(res_mat, BLUSPARAM)
        DataFrame(feature = names(ress),
                  cluster = clus,
                  sample_id = s)
    })
    if (length(sample_id) > 1L) {
        out <- do.call(rbind, out)
        out$cluster <- factor(out$cluster, levels = seq_len(max(as.integer(out$cluster))))
    } else out <- out[[1]]
    out
}


#' Find clusters on the Moran plot
#'
#' The Moran plot plots the value at each location on the x axis, and the
#' average of the neighbors of each locations on the y axis. Sometimes clusters
#' can be seen on the Moran plot, indicating different types of neighborhoods.
#'
#' @inheritParams bluster::clusterRows
#' @inheritParams calculateUnivariate
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
#' @importFrom SpatialFeatureExperiment localResults
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' library(bluster)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute moran plot
#' sfe <- runUnivariate(sfe, type = "moran.plot", features = rownames(sfe)[1],
#'                      exprs_values = "counts")
#' clusts <- clusterMoranPlot(sfe, rownames(sfe)[1],
#'                            BLUSPARAM = KmeansParam(2))
clusterMoranPlot <- function(sfe, features, BLUSPARAM, sample_id = NULL,
                             colGeometryName = NULL,
                             annotGeometryName = NULL, show_symbol = TRUE) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  use_col <- is.null(annotGeometryName) || !is.null(colGeometryName)
  if (is.null(colGeometryName) && is.null(annotGeometryName))
      features <- .symbol2id(sfe, features)
  out <- lapply(sample_id, function(s) {
    mps <- localResults(sfe, name = "moran.plot", features = features,
                        sample_id = s, colGeometryName = colGeometryName,
                        annotGeometryName = annotGeometryName)
    if (is.data.frame(mps)) o <- clusterRows(mps[,c("x", "wx")], BLUSPARAM)
    else
        o <- lapply(mps, function(mp) clusterRows(mp[,c("x", "wx")], BLUSPARAM))
    o <- as(o, "DataFrame")
    if (is.data.frame(mps)) names(o)[1] <- features
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
  if (show_symbol && any(features %in% rownames(sfe))) {
      ind <- features %in% rownames(sfe)
      features[ind] <- rowData(sfe)[featurs[ind], "symbol"]
  }
  out
}
