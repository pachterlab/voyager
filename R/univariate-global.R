# Univariate, from spdep

#' Calculate univariate global spatial autocorrelation
#'
#' Compute Moran's I or Geary's C on gene expression or numeric columns of
#' colData, colGeometry, or annotGeometry of a \code{SpatialFeatureExperiment}
#' object. Multithreading is supported when computing for numerous genes.
#'
#' @inheritParams spdep::moran
#' @param x For \code{calculateMoransI} and \code{calculateGearysC}, it can be a
#'   numeric matrix whose rows are features/genes, or a
#'   \code{SpatialFeatureExperiment} (SFE) object with such a matrix in an
#'   assay.
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param listw Weighted neighborhood graph as a \code{spdep} \code{listw}
#'   object.
#' @param features Genes (\code{calculate*} SFE method and \code{run*}) or
#'   numeric columns of \code{colData(x)} (\code{colData*}) or any
#'   \code{\link{colGeometry}} (\code{colGeometryM*}) or
#'   \code{\link{annotGeometry}} (\code{annotGeometry*}) for which the
#'   univariate metric is to be computed. Default to \code{NULL}. When
#'   \code{NULL}, then the metric is computed for all genes with the values in
#'   the assay specified in the argument \code{exprs_values}. This can be
#'   parallelized with the argument \code{BPPARAM}. For genes, if the column
#'   "symbol" is present in \code{rowData} and the row names of the SFE object
#'   are Ensembl IDs, then the gene symbol can be used and converted to IDs
#'   behind the scene. However, if one symbol matches multiple IDs, a warning
#'   will be given and the first match will be used.
#' @param exprs_values Integer scalar or string indicating which assay of x
#'   contains the expression values.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether
#'   and how computing the metric for numerous genes shall be parallelized.
#' @param name String specifying the name to be used to store the results in
#'   \code{rowData(x)}. If not already present in the name, then the
#'   \code{sample_id} will be appended to the name specified here separated by
#'   an underscore.
#' @param colGraphName Name of the listw graph in the SFE object that
#'   corresponds to entities represented by columns of the gene count matrix.
#'   Use \code{\link{colGraphNames}} to look up names of the available graphs
#'   for cells/spots. Note that for multiple \code{sample_id}s, it is assumed
#'   that all of them have a graph of this same name.
#' @param annotGraphName Name of the listw graph in the SFE object that
#'   corresponds to the \code{annotGeometry} of interest. Use
#'   \code{\link{annotGraphNames}} to look up names of available annotation
#'   graphs.
#' @param colGeometryName Name of a \code{colGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute the metric. Use
#'   \code{\link{colGeometryNames}} to look up names of the \code{sf} data
#'   frames associated with cells/spots.
#' @param annotGeometryName Name of a \code{annotGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute the metric. Use
#'   \code{\link{annotGeometryNames}} to look up names of the \code{sf} data
#'   frames associated with annotations.
#' @param sample_id Sample(s) in the SFE object whose cells/spots to use. Can be
#'   "all" to compute metric for all samples; the metric is computed separately
#'   for each sample.
#' @param ... Other arguments passed to \code{\link{moran}} or
#'   \code{\link{geary}}.
#' @return For \code{calculate*}, a \code{DataFrame} with two columns: The first
#'   one is I for Moran's I or C for Geary's C, and the second one is K for
#'   sample kurtosis. For the SFE method of \code{calculate*}, a third column
#'   indicating \code{sample_id} is added if more than one sample is indicated.
#'   For \code{run*}, a \code{SpatialFeatureExperiment} object with the Moran's
#'   I or Geary's C values added to a column of \code{rowData(x)}, whose name is
#'   specified in the \code{name} argument, with \code{sample_id} appended if
#'   applicable. For \code{colGeometry} and \code{annotGeometry}, the results
#'   are added to an attribute of the data frame called \code{featureData},
#'   which is a DataFrame analogous to \code{rowData} for the gene count matrix.
#'   New column names in \code{featureData} would follow the same rules as in
#'   \code{rowData}. For \code{colData}, the results can be accessed with the
#'   \code{colFeatureData} function.
#' @name calculateMoransI
#' @aliases calculateGearysC
#' @importFrom spdep moran geary Szero
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom SpatialFeatureExperiment SpatialFeatureExperiment
#' @importFrom SummarizedExperiment assay rowData<-
#' @importFrom SpatialFeatureExperiment colGraph annotGraph
#' @importFrom SingleCellExperiment colData rowData
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute Moran's I for vector or matrix
#' calculateMoransI(colData(sfe)$nCounts, listw = colGraph(sfe, "visium"))
#' # Add results to rowData, features are genes
#' sfe <- runMoransI(sfe, features = rownames(sfe)[1], exprs_values = "counts")
#' rowData(sfe)
#' # Specifically for colData
#' sfe <- colDataMoransI(sfe, "nCounts")
#' colFeatureData(sfe)
NULL

#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL) {
  .calc_univar_autocorr(x, listw, fun = moran, BPPARAM = BPPARAM,
                        n = length(listw$neighbours), S0 = Szero(listw),
                        zero.policy = zero.policy, returnDF = TRUE)
})

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL) {
  .calc_univar_autocorr(x, listw, fun = geary, BPPARAM = BPPARAM,
                        n = length(listw$neighbours),
                        n1 = length(listw$neighbours) - 1, S0 = Szero(listw),
                        zero.policy = zero.policy, returnDF = TRUE)
})

#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateMoransI, returnDF = TRUE))

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateGearysC, returnDF = TRUE))

#' @rdname calculateMoransI
#' @export
colDataMoransI <- .coldata_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
colDataGearysC <- .coldata_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

#' @rdname calculateMoransI
#' @export
colGeometryMoransI <- .colgeom_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
colGeometryGearysC <- .colgeom_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

#' @rdname calculateMoransI
#' @export
annotGeometryMoransI <- .annotgeom_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
annotGeometryGearysC <- .annotgeom_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

#' @rdname calculateMoransI
#' @export
runMoransI <- function(sfe, features, colGraphName = 1L, sample_id = NULL,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, name = "MoransI") {
  .sfe_univar_autocorr(sfe, features, colGraphName, sample_id, exprs_values,
                       calculateMoransI, BPPARAM, zero.policy, name = "MoransI")
}

#' @rdname calculateMoransI
#' @export
runGearysC <- function(sfe, features, colGraphName = 1L, sample_id = NULL,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, name = "MoransI") {
  .sfe_univar_autocorr(sfe, features, colGraphName, sample_id, exprs_values,
                       calculateGearysC, BPPARAM, zero.policy, name = "GearysC")
}

#' Permutation test for Moran's I and Geary's C
#'
#' Thin wrapper of \code{\link{moran.mc}} and \code{\link{geary.mc}} for easier
#' usage with the \code{SpatialFeatureExperiment} object and usage over multiple
#' genes. Multithreading is supported when computing for numerous genes.
#'
#' @inheritParams calculateMoransI
#' @inheritParams spdep::moran.mc
#' @param ... Other parameters passed to \code{\link{moran.mc}} or
#'   \code{\link{geary.mc}}.
#' @return For \code{calculateMoran/GearyMC}, a list of \code{mc.sim} objects.
#'   For the SFE method of \code{calculate*}, when more than one
#'   \code{sample_id} is specified, then a list of such lists, whose names are
#'   the \code{sample_id}s. For \code{runMoran/GearyMC}, the results are
#'   converted to a \code{DataFrame} and added to \code{rowData(x)}, and a SFE
#'   object with the added \code{rowData} is returned. For the
#'   colGeometry and annotGeometry versions, the results are added to the
#'   \code{featureData} attribute of the data frame of interest in a manner
#'   analogous to \code{rowData}. For \code{colData}, the results can be
#'   accessed with the \code{colFeatureData} function.
#' @importFrom spdep moran.mc geary.mc
#' @aliases calculateGearyMC
#' @name calculateMoranMC
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute Moran's I with Monte Carlo testing for vector or matrix
#' calculateMoranMC(colData(sfe)$nCounts, listw = colGraph(sfe, "visium"),
#'                  nsim = 100)
#' # Add results to rowData, features are genes
#' sfe <- runMoranMC(sfe, features = rownames(sfe)[1], exprs_values = "counts",
#'                   nsim = 100)
#' rowData(sfe)
#' # Specifically for colData
#' sfe <- colDataMoranMC(sfe, "nCounts", nsim = 100)
#' colFeatureData(sfe)
NULL

#' @rdname calculateMoranMC
#' @export
setMethod("calculateMoranMC", "ANY", function(x, listw, nsim,
                                              BPPARAM = SerialParam(),
                                              zero.policy = NULL,
                                              alternative = "greater", ...) {
  .calc_univar_autocorr(x, listw, fun = moran.mc, BPPARAM = BPPARAM,
                        nsim = nsim, zero.policy = zero.policy,
                        alternative = alternative, returnDF = FALSE, ...)
})

#' @rdname calculateMoranMC
#' @export
setMethod("calculateGearyMC", "ANY", function(x, listw, nsim,
                                              BPPARAM = SerialParam(),
                                              zero.policy = NULL,
                                              alternative = "greater", ...) {
  .calc_univar_autocorr(x, listw, fun = geary.mc, BPPARAM = BPPARAM,
                        nsim = nsim, zero.policy = zero.policy,
                        alternative = alternative, returnDF = FALSE, ...)
})

#' @rdname calculateMoranMC
#' @export
setMethod("calculateMoranMC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun_mc(calculateMoranMC))

#' @rdname calculateMoranMC
#' @export
setMethod("calculateGearyMC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun_mc(calculateGearyMC))

#' @rdname calculateMoranMC
#' @export
colDataMoranMC <- .coldata_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
colDataGearyMC <- .coldata_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

#' @rdname calculateMoranMC
#' @export
colGeometryMoranMC <- .colgeom_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
colGeometryGearyMC <- .colgeom_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

#' @rdname calculateMoranMC
#' @export
annotGeometryMoranMC <- .annotgeom_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
annotGeometryGearyMC <- .annotgeom_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

#' @rdname calculateMoranMC
#' @export
runMoranMC <- function(sfe, features, colGraphName = 1L, sample_id = NULL, nsim,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, alternative = "greater",
                       name = "MoranMC", ...) {
  .sfe_univar_mc(sfe, features, colGraphName, sample_id, nsim, exprs_values,
                 fun = calculateMoranMC, BPPARAM = BPPARAM,
                 zero.policy = zero.policy, alternative = alternative,
                 name = name, ...)
}

#' @rdname calculateMoranMC
#' @export
runGearyMC <- function(sfe, features, colGraphName = 1L, sample_id = NULL, nsim,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, alternative = "greater",
                       name = "GearyMC", ...) {
  .sfe_univar_mc(sfe, features, colGraphName, sample_id, nsim, exprs_values,
                 fun = calculateGearyMC, BPPARAM = BPPARAM,
                 zero.policy = zero.policy, alternative = alternative,
                 name = name, ...)
}

#' Spatial correlogram
#'
#' Still debating whether I should write the wrapper. It should be
#' straightforward to call sp.correlogram directly on single colData columns.
#' But for genes, there's more boilerplate. I suppose, for genes, it might be
#' cool to compute the correlogram for a bunch of genes and plot them in the
#' same plot, with the error bars, or cluster them. So I'll write the wrapper.
#'
#' @inheritParams calculateMoransI
#' @inheritParams spdep::sp.correlogram
#' @param ... Other arguments passed to \code{\link{sp.correlogram}}.
#' @importFrom spdep sp.correlogram
#' @return For \code{calculateCorrelogram}, a list of \code{spcor} objects, each
#'   element of which correslonds to a feature. or if multiple \code{sample_id}s
#'   are specified in the SFE method, a list of such lists whose names are the
#'   \code{sample_id}s. For \code{runCorrelogram}, the \code{res} field of the
#'   \code{spcor} is taken and put in a list column in \code{rowData(x)}, and
#'   the SFE object with the new \code{rowData} is returned. For the
#'   colGeometry and annotGeometry versions, the results are added to an
#'   attribute of the data frame of interest called \code{featureData}, in a
#'   manner analogous to \code{rowData}. For \code{colData}, the results can be
#'   accessed with the \code{colFeatureData} function.
#' @name calculateCorrelogram
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Compute correlogram for vector or matrix
#' calculateCorrelogram(colData(sfe)$nCounts, listw = colGraph(sfe, "visium"),
#'                      order = 5)
#' # Add results to rowData, features are genes
#' sfe <- runCorrelogram(sfe, features = rownames(sfe)[1], exprs_values = "counts",
#'                       order = 5)
#' rowData(sfe)
#' # Specifically for colData
#' sfe <- colDataCorrelogram(sfe, "nCounts", order = 5)
#' colFeatureData(sfe)
NULL

#' @rdname calculateCorrelogram
#' @export
setMethod("calculateCorrelogram", "ANY",
          function(x, listw, order = 1, method = "I", BPPARAM = SerialParam(),
                   zero.policy = NULL, ...) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  out <- bplapply(seq_len(nrow(x)), function(i) {
    sp.correlogram(listw$neighbours, var = x[i,], order = order, method = method,
                   zero.policy = zero.policy, ...)
  }, BPPARAM = BPPARAM)
  names(out) <- rownames(x)
  out
})

#' @rdname calculateCorrelogram
#' @export
setMethod("calculateCorrelogram", "SpatialFeatureExperiment",
          function(x, features, colGraphName = 1L, sample_id = NULL, order = 1,
                   method = "I", exprs_values = "logcounts",
                   BPPARAM = SerialParam(), zero.policy = NULL, ...) {
            .calc_univar_sfe_fun(calculateCorrelogram)(
              x, features, colGraphName, sample_id, exprs_values = exprs_values,
              BPPARAM = BPPARAM, zero.policy = zero.policy, order = order,
              method = method, ...)
          })

#' @rdname calculateCorrelogram
#' @export
colGeometryCorrelogram <- function(x, features, colGeometryName = 1L,
                                   colGraphName = 1L, sample_id = NULL,
                                   order = 1, method = "I",
                                   BPPARAM = SerialParam(), zero.policy = NULL,
                                   ...) {
  .colgeom_univar_fun(calculateCorrelogram, .correlogram2df,
                      name = paste("Correlogram", method, sep = "_"),
                      to_df_params = list(method = method))(
    x, features, colGeometryName, colGraphName, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
colDataCorrelogram <- function(x, features, colGraphName = 1L, sample_id = NULL,
                               order = 1, method = "I", BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .coldata_univar_fun(calculateCorrelogram, .correlogram2df,
                      name = paste("Correlogram", method, sep = "_"),
                      to_df_params = list(method = method))(
    x, features, colGraphName, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
annotGeometryCorrelogram <- function(x, features, annotGeometryName = 1L,
                                     annotGraphName = 1L, sample_id = NULL,
                                     order = 1, method = "I",
                                     BPPARAM = SerialParam(),
                                     zero.policy = NULL, ...) {
  .annotgeom_univar_fun(calculateCorrelogram, .correlogram2df,
                        name = paste("Correlogram", method, sep = "_"),
                        to_df_params = list(method = method))(
    x, features, annotGeometryName, annotGraphName, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
runCorrelogram <- function(sfe, features, colGraphName = 1L, sample_id = NULL,
                           order = 1, method = "I", exprs_values = "logcounts",
                           BPPARAM = SerialParam(), zero.policy = NULL,
                           name = paste("Correlogram", method, sep = "_"), ...) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  for (s in sample_id) {
    out <- calculateCorrelogram(sfe, features, colGraphName, s, order,
                                method, exprs_values, BPPARAM, zero.policy, ...)
    out <- .correlogram2df(out, name, method)
    features <- .symbol2id(sfe, features)
    rownames(out) <- features
    out <- .add_name_sample_id(out, s)
    rowData(sfe)[features, names(out)] <- out
  }
  sfe
}

#' Find clusters of correlogram patterns
#'
#' Cluster the correlograms to find patterns in length scales of spatial
#' autocorrelation. All the correlograms clustered must be computed with the
#' same method and have the same number of lags.
#'
#' @inheritParams clusterMoranPlot
#' @inheritParams calculateCorrelogram
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
#' sfe <- runCorrelogram(sfe, features = rownames(sfe)[inds],
#'                       exprs_values = "counts", order = 5)
#' clust <- clusterCorrelograms(sfe, features = rownames(sfe)[inds],
#'                              BLUSPARAM = KmeansParam(2))
clusterCorrelograms <- function(sfe, features, BLUSPARAM, sample_id = NULL,
                                method = "I",
                                name = paste("Correlogram", method, sep = "_"),
                                colGeometryName = NULL,
                                annotGeometryName = NULL, show_symbol = TRUE) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
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
