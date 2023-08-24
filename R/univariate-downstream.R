#' Find clusters of correlogram patterns
#'
#' Cluster the correlograms to find patterns in length scales of spatial
#' autocorrelation. All the correlograms clustered must be computed with the
#' same method and have the same number of lags. Correlograms are clustered
#' jointly across samples.
#'
#' @inheritParams clusterMoranPlot
#' @inheritParams plotCorrelogram
#' @inheritParams plotCorrelogram
#' @inheritParams plotDimLoadings
#' @param sfe A \code{SpatialFeatureExperiment} object with correlograms
#'   computed for features of interest.
#' @param features Features whose correlograms to cluster.
#' @return A data frame with 3 columns: \code{feature} for the features,
#'   \code{cluster} a factor for cluster membership of the features within each
#'   sample, and \code{sample_id} for the sample.
#' @concept Downstream analyses of univariate spatial results
#' @export
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SFEData)
#' library(bluster)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' inds <- c(1, 3, 4, 5)
#' sfe <- runUnivariate(sfe,
#'     type = "sp.correlogram",
#'     features = rownames(sfe)[inds],
#'     exprs_values = "counts", order = 5
#' )
#' clust <- clusterCorrelograms(sfe,
#'     features = rownames(sfe)[inds],
#'     BLUSPARAM = KmeansParam(2)
#' )
clusterCorrelograms <- function(sfe, features, BLUSPARAM, sample_id = "all",
                                method = "I", colGeometryName = NULL,
                                annotGeometryName = NULL, reducedDimName = NULL,
                                show_symbol = deprecated(),
                                swap_rownames = NULL, name = "sp.correlogram") {
    l <- .deprecate_show_symbol("clusterCorrelograms", show_symbol, swap_rownames)
    show_symbol <- l[[1]]; swap_rownames <- l[[2]]

    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    name <- paste(name, method, sep = "_")
    res <- lapply(sample_id, function(s) {
        ress <- .get_feature_metadata(sfe, features, name, s, colGeometryName,
            annotGeometryName, reducedDimName,
            show_symbol = show_symbol, swap_rownames = swap_rownames
        )
        if (method %in% c("I", "C")) {
            # First column is the metric, second column expectation, third is variance
            ress <- lapply(ress, function(r) r[, 1])
        }
        res_mat <- do.call(rbind, ress)
        rownames(res_mat) <- paste(names(ress), s, sep = "__")
        res_mat
    })
    res <- do.call(rbind, res)
    nn <- do.call(rbind, strsplit(rownames(res), "__"))
    clus <- clusterRows(res, BLUSPARAM)
    out <- data.frame(
        feature = nn[,1],
        cluster = clus,
        sample_id = nn[,2]
    )
    rownames(out) <- NULL
    out$cluster <- factor(out$cluster,
                          levels = seq_len(max(as.integer(out$cluster))))
    out
}

#' Find clusters on the Moran plot
#'
#' The Moran plot plots the value at each location on the x axis, and the
#' average of the neighbors of each locations on the y axis. Sometimes clusters
#' can be seen on the Moran plot, indicating different types of neighborhoods.
#'
#' @inheritParams bluster::clusterRows
#' @inheritParams plotCorrelogram
#' @inheritParams plotDimLoadings
#' @param sfe A \code{SpatialFeatureExperiment} object with Moran plot computed
#'   for the feature of interest. If the Moran plot for that feature has not
#'   been computed for that feature in this sample_id, it will be calculated and
#'   stored in \code{rowData}. See \code{\link{calculateUnivariate}}.
#' @param colGeometryName Name of colGeometry from which to look for features.
#' @param annotGeometryName Name of annotGeometry from which to look for
#'   features.
#' @param features Features whose Moran plot are to be cluster. Features whose
#'   Moran plots have not been computed will be skipped, with a warning.
#' @return A data frame each column of which is a factor for cluster
#'   membership of each feature. The column names are the features.
#' @importFrom bluster clusterRows
#' @importFrom methods as
#' @importFrom SpatialFeatureExperiment localResults
#' @export
#' @concept Downstream analyses of univariate spatial results
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' library(bluster)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute moran plot
#' sfe <- runUnivariate(sfe,
#'     type = "moran.plot", features = rownames(sfe)[1],
#'     exprs_values = "counts"
#' )
#' clusts <- clusterMoranPlot(sfe, rownames(sfe)[1],
#'     BLUSPARAM = KmeansParam(2)
#' )
clusterMoranPlot <- function(sfe, features, BLUSPARAM, sample_id = "all",
                             colGeometryName = NULL,
                             annotGeometryName = NULL,
                             show_symbol = deprecated(), swap_rownames = NULL) {
    l <- .deprecate_show_symbol("clusterMoranPlot", show_symbol, swap_rownames)
    show_symbol <- l[[1]]; swap_rownames <- l[[2]]

    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    use_col <- is.null(annotGeometryName) || !is.null(colGeometryName)
    if (is.null(colGeometryName) && is.null(annotGeometryName)) {
        features <- .symbol2id(sfe, features, swap_rownames)
    }
    out <- lapply(sample_id, function(s) {
        mps <- localResults(sfe,
            name = "moran.plot", features = features,
            sample_id = s, colGeometryName = colGeometryName,
            annotGeometryName = annotGeometryName
        )
        if (is.data.frame(mps)) {
            o <- clusterRows(mps[, c("x", "wx")], BLUSPARAM)
        } else {
            o <- lapply(mps, function(mp)
                clusterRows(mp[, c("x", "wx")], BLUSPARAM))
        }
        o <- as.data.frame(o)
        if (is.data.frame(mps)) names(o)[1] <- features
        o$sample_id <- s
        if (use_col) {
            row.names(o) <- colnames(sfe)[colData(sfe)$sample_id == s]
        }
        # What if some features don't have the Moran Plot computed
        features_absent <- setdiff(features, names(o))
        if (length(features_absent) && length(sample_id) > 1L) {
            for (f in features_absent) {
                o[[f]] <- NA
            }
            o <- o[, c("sample_id", features)]
        }
        o
    })
    if (length(sample_id) > 1L) {
        out <- do.call(rbind, out)
        for (f in features) {
            if (all(is.na(out[[f]]))) {
                out[[f]] <- NULL
            } else {
                out[[f]] <- factor(out[[f]],
                                   levels = seq_len(max(as.integer(out[[f]]))))
            }
        }
    } else {
        out <- out[[1]]
    }
    if (show_symbol && any(features %in% rownames(sfe)) &&
        "symbol" %in% names(rowData(sfe))) {
        ind <- features %in% rownames(sfe)
        features[ind] <- rowData(sfe)[features[ind], "symbol"]
    }
    out
}

#' Cluster variograms of multiple features
#'
#' This function clusters variograms of features across samples to find patterns
#' in decays in spatial autocorrelation. The fitted variograms are clustered as
#' different samples can have different distance bins.
#'
#' @inheritParams clusterCorrelograms
#' @param n Number of points on the fitted variogram line.
#' @return A data frame with 3 columns: \code{feature} for the features,
#'   \code{cluster} a factor for cluster membership of the features within each
#'   sample, and \code{sample_id} for the sample.
#' @export
#' @concept Downstream analyses of univariate spatial results
#' @examples
#' library(SFEData)
#' library(scater)
#' library(bluster)
#' library(Matrix)
#' sfe <- McKellarMuscleData()
#' sfe <- logNormCounts(sfe)
#' # Just the highly expressed genes
#' gs <- order(Matrix::rowSums(counts(sfe)), decreasing = TRUE)[1:10]
#' genes <- rownames(sfe)[gs]
#'
#' sfe <- runUnivariate(sfe, "variogram", features = genes)
#' clusts <- clusterVariograms(sfe, genes, BLUSPARAM = HclustParam(),
#' swap_rownames = "symbol")
#'
#' # Plot the clustering
#' plotVariogram(sfe, genes, color_by = clusts, group = "feature",
#' use_lty = FALSE, swap_rownames = "symbol", show_np = FALSE)
#'
clusterVariograms <- function(sfe, features, BLUSPARAM, n = 20,
                              sample_id = "all", colGeometryName = NULL,
                              annotGeometryName = NULL, reducedDimName = NULL,
                              swap_rownames = NULL, name = "variogram") {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    rlang::check_installed("gstat")
    show_symbol <- !is.null(swap_rownames)
    ress <- lapply(sample_id, function(s) {
        .get_feature_metadata(sfe, features, name, s, colGeometryName,
                              annotGeometryName, reducedDimName,
                              show_symbol = show_symbol,
                              swap_rownames = swap_rownames
        )
    })
    max_dists <- lapply(ress, function(res)
        vapply(res, function(r) max(r$exp_var$dist),
               FUN.VALUE = numeric(1))
    )
    max_dists <- unlist(max_dists)
    dist_use <- min(max_dists)

    res_mat <- lapply(seq_along(ress), function(i) {
        var_lines <- lapply(seq_along(ress[[i]]), function(j) {
            m <- ress[[i]][[j]]$var_model
            l <- gstat::variogramLine(m, n = n, maxdist = dist_use)
            l$gamma
        })
        m <- do.call(rbind, var_lines)
        rownames(m) <- paste(names(ress[[i]]), sample_id[i], sep = "__")
        m
    })
    res_mat <- do.call(rbind, res_mat)
    nn <- do.call(rbind, strsplit(rownames(res_mat), "__"))
    clus <- clusterRows(res_mat, BLUSPARAM)
    out <- data.frame(
        feature = nn[,1],
        cluster = clus,
        sample_id = nn[,2]
    )
    rownames(out) <- NULL
    out$cluster <- factor(out$cluster,
                          levels = seq_len(max(as.integer(out$cluster))))
    out
}
