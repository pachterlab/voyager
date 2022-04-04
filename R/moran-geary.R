# Univariate, from spdep

#' Calculate univariate spatial autocorrelation
#'
#' Compute Moran's I or Geary's C on gene expression or numeric columns of
#' colData, colGeometry, or annotGeometry of a \code{SpatialFeatureExperiment}
#' object. Multithreading is supported when computing for numerous genes.
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
#'   and how computing Moran's I or Geary's C for numerous genes shall be
#'   parallelized.
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

.calc_univar_autocorr <- function(x, listw, fun, BPPARAM, returnDF = TRUE, ...) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  out_list <- bplapply(seq_len(nrow(x)), function(i) {
    fun(x[i,], listw, ...)
  }, BPPARAM = BPPARAM)
  if (returnDF) {
    out_list <- lapply(out_list, unlist, use.names = TRUE)
    out <- Reduce(rbind, out_list)
    rownames(out) <- names(out_list)
    out <- DataFrame(out)
  }
  return(out)
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

.calc_univar_sfe_fun <- function(fun) {
  function(x, colGraphName, features, sample_id, exprs_values = "logcounts",
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    # Am I sure that I want to use logcounts as the default?
    if (!all(features %in% rownames(x))) {
      features <- intersect(features, rownames(x))
      if (!length(features)) {
        stop("None of the specified genes/features are found in the SFE object.")
      }
    }
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    mat <- assay(x, exprs_values)[features, colData(x)$sample_id %in% sample_id]
    fun(mat, listw_use, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
  }
}
#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateMoransI))

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateGearysC))

.df_univar_autocorr <- function(df, listw, features, fun, BPPARAM,
                                zero.policy, ...) {
  mat <- t(as.matrix(df[, features]))
  if (anyNA(mat)) {
    stop("Only numeric columns without NA (within the sample_id) can be used.")
  }
  fun(mat, listw, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
}

.coldata_univar_fun <- function(fun) {
  function(x, colGraphName, features, sample_id, BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    .df_univar_autocorr(colData(x)[colData(x)$sample_id %in% sample_id,],
                        listw_use, features, fun,
                        BPPARAM, zero.policy, ...)
  }
}

#' @rdname calculateMoransI
#' @export
colDataMoransI <- .coldata_univar_fun(calculateMoransI)

#' @rdname calculateMoransI
#' @export
colDataGearysC <- .coldata_univar_fun(calculateGearysC)

.colgeom_univar_fun <- function(fun) {
  function(x, colGeometryName, colGraphName, features, sample_id,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    .df_univar_autocorr(colGeometry(x, type = colGeometryName,
                                    sample_id = sample_id),
                        listw_use, features, fun, BPPARAM,
                        zero.policy, ...)
  }
}

#' @rdname calculateMoransI
#' @export
colGeometryMoransI <- .colgeom_univar_fun(calculateMoransI)

#' @rdname calculateMoransI
#' @export
colGeometryGearysC <- .colgeom_univar_fun(calculateGearysC)

.annotgeom_univar_fun <- function(fun) {
  function(x, annotGeometryName, annotGraphName, features, sample_id,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    listw_use <- annotGraph(x, type = annotGraphName, sample_id = sample_id)
    .df_univar_autocorr(annotGeometry(x, type = annotGeometryName,
                                      sample_id = sample_id),
                        listw_use, features, fun, BPPARAM,
                        zero.policy, ...)
  }
}

#' @rdname calculateMoransI
#' @export
annotGeometryMoransI <- .annotgeom_univar_fun(calculateMoransI)

#' @rdname calculateMoransI
#' @export
annotGeometryGearysC <- .annotgeom_univar_fun(calculateGearysC)

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

#' Permutation test for Moran's I and Geary's C
#'
#' Thin wrapper of \code{\link{moran.mc}} and \code{\link{geary.mc}} for easier
#' usage with the \code{SpatialFeatureExperiment} object and usage over multiple
#' genes. Multithreading is supported when computing for numerous genes.
#'
#' @inheritParams calculateMoransI
#' @inheritParams spdep::moran.mc
#' @param ... Other parameters passed to \code{\link{moran.mc}} or
#' \code{\link{geary.mc}}.
#' @return For \code{calculateMoran/GearyMC}, a list of \code{mc.sim} objects.
#' For \code{runMoran/GearyMC}, the results are converted to a \code{DataFrame}
#' and added to \code{rowData(x)}, and a SFE object with the added \code{rowData}
#' is returned.
#' @importFrom spdep moran.mc geary.mc
#' @aliases calculateGearyMC
#' @name calculateMoranMC
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

# For the `nsim` and `alternative` arguments
.calc_univar_sfe_fun_mc <- function(fun) {
  function(x, colGraphName, features, sample_id, nsim, exprs_values = "logcounts",
           BPPARAM = SerialParam(), zero.policy = NULL,
           alternative = "greater", ...) {
    .calc_univar_sfe_fun(fun)(x, colGraphName, features, sample_id,
                              exprs_values, BPPARAM, zero.policy, nsim = nsim,
                              alternative = alternative, ...)
  }
}

#' @rdname calculateMoranMC
#' @export
setMethod("calculateMoranMC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun_mc(calculateMoranMC))

#' @rdname calculateMoranMC
#' @export
setMethod("calculateGearyMC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun_mc(calculateGearyMC))

.coldata_univar_fun_mc <- function(fun) {
  function(x, colGeometryName, colGraphName, features, sample_id, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .coldata_univar_fun(fun)(x, colGeometryName, colGraphName, features,
                             sample_id, BPPARAM, zero.policy, nsim = nsim,
                             alternative = alternative, ...)
  }
}

#' @rdname calculateMoranMC
#' @export
colDataMoranMC <- .coldata_univar_fun_mc(calculateMoranMC)

#' @rdname calculateMoranMC
#' @export
colDataGearyMC <- .coldata_univar_fun_mc(calculateGearyMC)

.colgeom_univar_fun_mc <- function(fun) {
  function(x, colGeometryName, colGraphName, features, sample_id, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .colgeom_univar_fun(fun)(x, colGeometryName, colGraphName, features,
                             sample_id, BPPARAM, zero.policy, nsim = nsim,
                             alternative = alternative, ...)
  }
}
#' @rdname calculateMoranMC
#' @export
colGeometryMoranMC <- .colgeom_univar_fun_mc(calculateMoranMC)

#' @rdname calculateMoranMC
#' @export
colGeometryGearyMC <- .colgeom_univar_fun_mc(calculateGearyMC)

.annotgeom_univar_fun_mc <- function(fun) {
  function(x, annotGeometryName, annotGraphName, features, sample_id, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .annotgeom_univar_fun(fun)(x, annotGeometryName, annotGraphName, features,
                               sample_id, BPPARAM, zero.policy, nsim = nsim,
                               alternative = alternative, ...)
  }
}

#' @rdname calculateMoranMC
#' @export
annotGeometryMoranMC <- .annotgeom_univar_fun_mc(calculateMoranMC)

#' @rdname calculateMoranMC
#' @export
annotGeometryGearyMC <- .annotgeom_univar_fun_mc(calculateGearyMC)

.sfe_univar_mc <- function(x, colGraphName, features, sample_id, nsim,
                           exprs_values, fun, BPPARAM, zero.policy, alternative,
                           name, ...) {
  out <- fun(x, colGraphName, features, sample_id, nsim, exprs_values, BPPARAM,
             zero.policy, alternative, ...)
  # Convert results to DFrame. I'll write my own plotting function based on ggplot2.
  out <- lapply(out, function(o) {
    o$res <- I(list(o$res))
    DataFrame(unclass(o))
  })
  rns <- names(out)
  out <- Reduce(out, rbind)
  rownames(out) <- rns
  colnames(out) <- paste(name, colnames(out), sep = "_")
  if (length(sampleIDs(x)) > 1L) {
    names(out) <- paste(names(out), sample_id, sep = "_")
  }
  rowData(x)[features, names(out)] <- out
  x
}

#' @rdname calculateMoranMC
#' @export
runMoranMC <- function(x, colGraphName, features, sample_id, nsim,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, alternative = "greater",
                       name = "MoranMC", ...) {
  .sfe_univar_mc(x, colGraphName, features, sample_id, nsim, exprs_values,
                 fun = calculateMoranMC, BPPARAM = BPPARAM,
                 zero.policy = zero.policy, alternative = alternative,
                 name = name, ...)
}

#' @rdname calculateMoranMC
#' @export
runGearyMC <- function(x, colGraphName, features, sample_id, nsim,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, alternative = "greater",
                       name = "GearyMC", ...) {
  .sfe_univar_mc(x, colGraphName, features, sample_id, nsim, exprs_values,
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
#' @return For \code{calculateCorrelogram} and the \code{colData},
#'   \code{colGeometry}, and \code{annotGeometry} versions, a list of
#'   \code{spcor} objects, each element of which correslonds to a feature. For
#'   \code{runCorrelogram}, the \code{res} field of the \code{spcor} is taken
#'   and put in a list column in \code{rowData(x)}, and the SFE object with the
#'   new \code{rowData} is returned.
#' @name calculateCorrelogram
NULL

#' @rdname calculateCorrelogram
#' @export
setMethod("calculateCorrelogram", "ANY",
          function(x, listw, order = 1, method = "I", BPPARAM = SerialParam(),
                   zero.policy = NULL, ...) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  bplapply(seq_len(nrow(x)), function(i) {
    sp.correlogram(listw$neighbours, var = x[i,], order = order, method = method,
                   zero.policy = zero.policy, ...)
  }, BPPARAM = BPPARAM)
})

#' @rdname calculateCorrelogram
#' @export
setMethod("calculateCorrelogram", "SpatialFeatureExperiment",
          function(x, colGraphName, features, sample_id, order = 1,
                   method = "I", exprs_values = "logcounts",
                   BPPARAM = SerialParam(), zero.policy = NULL, ...) {
            .calc_univar_sfe_fun(calculateCorrelogram)(
              x, colGraphName, features, sample_id, exprs_values = exprs_values,
              BPPARAM = BPPARAM, zero.policy = zero.policy, order = order,
              method = method, ...)
          })

#' @rdname calculateCorrelogram
#' @export
colGeometryCorrelogram <- function(x, colGeometryName, colGraphName, features,
                               sample_id, order = 1, method = "I",
                               BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .colgeom_univar_fun(calculateCorrelogram)(
    x, colGeometryName, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
colDataCorrelogram <- function(x, colGraphName, features, sample_id,
                               order = 1, method = "I", BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .coldata_univar_fun(calculateCorrelogram)(
    x, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
annotGeometryCorrelogram <- function(x, annotGeometryName, annotGraphName,
                                     features, sample_id, order = 1,
                                     method = "I", BPPARAM = SerialParam(),
                                     zero.policy = NULL, ...) {
  .annotgeom_univar_fun(calculateCorrelogram)(
    x, annotGeometryName, annotGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
runCorrelogram <- function(x, colGraphName, features, sample_id, order = 1,
                           method = "I", exprs_values = "logcounts",
                           BPPARAM = SerialParam(), zero.policy = NULL,
                           name = paste("Correlogram", method, sep = "_"), ...) {
  out <- calculateCorrelogram(x, colGraphName, features, sample_id, order,
                              method, exprs_values, BPPARAM, zero.policy, ...)
  out <- lapply(out, function(o) I(list(o$res)))
  out_df <- DataFrame(res = out, row.names = names(out))
  names(out_df) <- name
  if (length(sampleIDs(x)) > 1L) {
    names(out_df) <- paste(names(out_df), sample_id, sep = "_")
  }
  rowData(x)[features, names(out_df)] <- out_df
  x
}
