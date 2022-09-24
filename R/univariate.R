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
#' @param ... Other arguments passed to S4 method (for convenience wrappers like
#'   \code{calculateMoransI}) or method used to compute metrics as specified by
#'   the argument \code{type} (as in more general functions like
#'   \code{calculateUnivariate}). See documentation in the \code{spdep} package
#'   for the latter.
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
#' @name calculateUnivariate
#' @aliases calculateMoransI
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

#' @rdname calculateUnivariate
#' @export
setMethod("calculateUnivariate", "ANY",
          function(x, listw, type = c("moran", "geary", "moran.mc", "geary.mc",
                                      "moran.test", "geary.test", "globalG.test",
                                      "sp.correlogram", "moran.plot", "localmoran",
                                      "localmoran_perm", "localC", "localC_perm",
                                      "loclaG", "localG_perm", "LOSH", "LOSH.mc",
                                      "gwss"),
                   BPPARAM = SerialParam(),
                   zero.policy = NULL, ...) {
              type <- match.arg(type)
              fun <- match.fun(type)
              obscure_args <- switch(type,
                                     moran = c("n", "S0"),
                                     geary = c("n", "n1", "S0"))
              defaults <- .obscure_arg_defaults(listw, type)
              other_args <- list(...)
              defaults_use <- defaults[setdiff(names(defaults), other_args)]
              all_args <- list(x = x, listw = listw, fun = fun,
                               BPPARAM = BPPARAM, returnDF = returnDF,
                               zero.policy = zero.policy)
              all_args <- c(all_args, other_args, defaults_use)
              out <- do.call(.calc_univar, all_args)
              out
          })

#' @rdname calculateUnivariate
#' @export
setMethod("calculateUnivariate", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun())

#' @rdname calculateUnivariate
#' @export
setMethod("calculateMoransI", "ANY",
          function(x, ..., BPPARAM = SerialParam(), zero.policy = NULL)
              calculateUnivariate(x, type = "moran", BPPARAM = BPPARAM,
                                  zero.policy = zero.policy, ...)
          )

#' @rdname calculateUnivariate
#' @export
setMethod("calculateMoransI", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(type = "moran"))

#' @rdname calculateUnivariate
#' @export
colDataUnivariate <- .coldata_univar_fun()

#' @rdname calculateUnivariate
#' @export
colDataMoransI <- .coldata_univar_fun(type = "moran")

#' @rdname calculateUnivariate
#' @export
colGeometryUnivariate <- .colgeom_univar_fun()

#' @rdname calculateUnivariate
#' @export

colGeometryMoransI <- .colgeom_univar_fun(type = "moran")

#' @rdname calculateUnivariate
#' @export
annotGeometryUnivariate <- .annotgeom_univar_fun()

#' @rdname calculateUnivariate
#' @export
annotGeometryMoransI <- .annotgeom_univar_fun(type = "moran")

#' @rdname calculateUnivariate
#' @export
runUnivariate <- .sfe_univar_fun()

#' @rdname calculateUnivariate
#' @export
runMoransI <- .sfe_univar_fun(type = "moran")
