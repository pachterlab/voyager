# Univariate, from spdep
# 2. For moran.mc, results can be stored in rowData as well
# 3. moran.plot: The output dataframe can be stored in colData. Also options to
# just return the results. Just like calculatePCA vs runPCA in scater
# 4. localmoran: Output dataframe can be stored in colData
# 5. Clustering based on moran.plot results
# 6. Find high-high, low-low, high-low, low-high regions from moran.plot and localmoran
# 7. Getis-Ord Gi and Gi*
# 8. My own plotting function for moran.plot, with ggplot2
# 9. Plotting with divergent palette
# 10. What to do with the image when using geom_sf
# 11. LOSH
# 12. Correlogram (sp.correlogram). Results can be stored in rowData or separately

#' Calculate univariate spatial autocorrelation
#'
#' Compute Moran's I or Geary's C on gene expression or numeric columns of
#' colData, colGeometry, or annotGeometry of a \code{SpatialFeatureExperiment}
#' object.
#'
#' @inheritParams spdep::moran
#' @param x For \code{calculateMoransI} and \code{calculateGearysC}, it can be a
#'   numeric matrix whose rows are features/genes, or a
#'   \code{SpatialFeatureExperiment} (SFE) object with such a matrix in an
#'   assay. For \code{runMoransI} and \code{runGearysC}, and the \code{colData},
#'   \code{colGeometry}, and \code{annotGeometry} versions, it must be a
#'   \code{SpatialFeatureExperiment} object.
#' @param listw Weighted neighborhood graph as a \code{spdep} \code{listw}
#'   object.
#' @param features Genes (\code{calculateMoransI/GearysC} SFE method and
#'   \code{runMoransI/GearysC}) or numeric columns of \code{colData(x)}
#'   (\code{colDataMoransI/GearysC}) or any \code{\link{colGeometry}}
#'   (\code{colGeometryMoransI/GearysC}) or \code{\link{annotGeometry}}
#'   (\code{annotGeometryMoransI/GearysC}) for which Moran's I or Geary's C is
#'   to be computed. Default to \code{NULL}. When \code{NULL}, then Moran's I or
#'   Geary's C is computed for all genes with the values in the assay specified
#'   in the argument \code{exprs_values}. This can be parallelized with the
#'   argument \code{BPPARAM}.
#' @param exprs_values Integer scalar or string indicating which assay of x
#'   contains the expression values.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether
#'   and how computing Moran's I for numerous genes shall be parallelized.
#' @param name String specifying the name to be used to store the results in
#'   \code{rowData(x)} for \code{runMoransI} and \code{runGearysC}. If the SFE
#'   object has more than one \code{sample_id}, then the \code{sample_id} will
#'   be appended to the name specified here separated by an underscore.
#' @param colGraphName Name of the listw graph in the SFE object that
#'   corresponds to entities represented by columns of the gene count matrix.
#'   Use \code{\link{colGraphNames}} to look up names of the available graphs
#'   for cells/spots.
#' @param annotGraphName Name of the listw graph in the SFE object that
#'   corresponds to the \code{annotGeometry} of interest. Use
#'   \code{\link{annotGraphNames}} to look up names of available annotation
#'   graphs.
#' @param colGeometryName Name of a \code{colGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute Moran's I or
#'   Geary's C. Use \code{\link{colGeometryNames}} to look up names of the
#'   \code{sf} data frames associated with cells/spots.
#' @param annotGeometryName Name of a \code{annotGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute Moran's I or
#'   Geary's C. Use \code{\link{annotGeometryNames}} to look up names of the
#'   \code{sf} data frames associated with annotations.
#' @param sample_id Sample in the SFE object whose cells/spots to use to
#'   calculate Moran's I or Geary's C.
#' @return For \code{calculate*} and the \code{colData}, \code{colGeometry}, and
#'   \code{annotGeometry} versions, a \code{DataFrame} with two columns: The
#'   first one is I for Moran's I or C for Geary's C, and the second one is K
#'   for sample kurtosis. For \code{run*}, a \code{SpatialFeatureExperiment}
#'   object with the Moran's I or Geary's C values added to a column of
#'   \code{rowData(x)}, whose name is specified in the \code{name} argument,
#'   with \code{sample_id} appended if applicable.
#' @name calculateMoransI
#' @aliases calculateGearysC
#' @importFrom spdep moran geary
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom SpatialFeatureExperiment SpatialFeatureExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialFeatureExperiment colGraph annotGraph
#' @importFrom SingleCellExperiment colData
NULL

.calc_univar_autocorr <- function(x, listw, fun, BPPARAM, ...) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  out_list <- bplapply(seq_len(nrow(x)), function(i) {
    unlist(fun(x[i,], listw, ...), use.names = TRUE)
  }, BPPARAM = BPPARAM)
  out <- Reduce(rbind, out_list)
  rownames(out) <- names(out_list)
  DataFrame(out)
}

#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL) {
  .calc_univar_autocorr(x, listw, fun = moran, BPPARAM = BPPARAM,
                        n = length(listw$neighbours), S0 = Szero(listw),
                        zero.policy = zero.policy)
})

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL) {
  .calc_univar_autocorr(x, listw, fun = geary, BPPARAM = BPPARAM,
                        n = length(listw$neighbours),
                        n1 = length(listw$neighbours) - 1, S0 = Szero(listw),
                        zero.policy = zero.policy)
})

#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "SpatialFeatureExperiment",
          function(x, colGraphName, features, sample_id,
                   exprs_values = "logcounts", BPPARAM = SerialParam(),
                   zero.policy = NULL) {
            # Am I sure that I want to use logcounts as the default?
            listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
            mat <- assay(x, exprs_values)[,colData(x)$sample_id %in% sample_id]
            calculateMoransI(mat, listw_use, BPPARAM, zero.policy)
          })

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "SpatialFeatureExperiment",
          function(x, colGraphName, features, sample_id,
                   exprs_values = "logcounts", BPPARAM = SerialParam(),
                   zero.policy = NULL) {
            listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
            mat <- assay(x, exprs_values)[,colData(x)$sample_id %in% sample_id]
            calculateGearysC(mat, listw_use, BPPARAM, zero.policy)
          })

.df_univar_autocorr <- function(df, listw, features, fun, BPPARAM, zero.policy) {
  mat <- t(as.matrix(df[, features]))
  if (anyNA(mat)) {
    stop("Only numeric columns without NA (within the sample_id) can be used.")
  }
  fun(mat, listw, BPPARAM, zero.policy)
}

#' @rdname calculateMoransI
#' @export
colDataMoransI <- function(x, colGraphName, features, sample_id,
                           BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
  .df_univar_autocorr(colData(x)[colData(x)$sample_id %in% sample_id],
                      listw_use, features, calculateMoransI,
                      BPPARAM, zero.policy)
}

#' @rdname calculateMoransI
#' @export
colDataGearysC <- function(x, colGraphName, features, sample_id,
                           BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
  .df_univar_autocorr(colData(x)[colData(x)$sample_id %in% sample_id],
                      listw_use, features, calculateGearysC,
                      BPPARAM, zero.policy)
}

#' @rdname calculateMoransI
#' @export
colGeometryMoransI <- function(x, colGeometryName, colGraphName, features,
                               sample_id,
                               BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
  .df_univar_autocorr(colGeometry(x, type = colGeometryName,
                                  sample_id = sample_id),
                      listw_use, features, calculateMoransI, BPPARAM,
                      zero.policy)
}

#' @rdname calculateMoransI
#' @export
colGeometryGearysC <- function(x, colGeometryName, colGraphName, features,
                               sample_id,
                               BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
  .df_univar_autocorr(colGeometry(x, type = colGeometryName,
                                  sample_id = sample_id),
                      listw_use, features, calculateGearysC, BPPARAM,
                      zero.policy)
}

#' @rdname calculateMoransI
#' @export
annotGeometryMoransI <- function(x, annotGeometryName, annotGraphName, features,
                                 sample_id,
                                 BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- annotGraph(x, type = annotGraphName, sample_id = sample_id)
  .df_univar_autocorr(annotGeometry(x, type = annotGeometryName,
                                    sample_id = sample_id),
                      listw_use, features, calculateMoransI, BPPARAM,
                      zero.policy)
}

#' @rdname calculateMoransI
#' @export
annotGeometryGearysC <- function(x, annotGeometryName, annotGraphName, features,
                                 sample_id,
                                 BPPARAM = SerialParam(), zero.policy = NULL) {
  listw_use <- annotGraph(x, type = annotGraphName, sample_id = sample_id)
  .df_univar_autocorr(annotGeometry(x, type = annotGeometryName,
                                    sample_id = sample_id),
                      listw_use, features, calculateGearysC, BPPARAM,
                      zero.policy)
}

.sfe_univar_autocorr <- function(x, colGraphName, features, sample_id,
                                 exprs_values, fun, BPPARAM, zero.policy,
                                 name) {
  out <- fun(x, colGraphName, features, sample_id, exprs_values, BPPARAM,
             zero.policy)
  names(out)[1] <- name
  if (length(sampleIDs(x)) > 1L) {
    names(out) <- paste(names(out), sample_id, sep = "_")
  }
  rowData(x)[features, names(out)] <- out
  x
}

#' @rdname calculateMoransI
#' @export
runMoransI <- function(x, colGraphName, features, sample_id,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, name = "MoransI") {
  .sfe_univar_autocorr(x, colGraphName, features, sample_id, exprs_values,
                       calculateMoransI, BPPARAM, zero.policy, name = "MoransI")
}

#' @rdname calculateMoransI
#' @export
runGearysC <- function(x, colGraphName, features, sample_id,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, name = "MoransI") {
  .sfe_univar_autocorr(x, colGraphName, features, sample_id, exprs_values,
                       calculateGearysC, BPPARAM, zero.policy, name = "GearysC")
}
