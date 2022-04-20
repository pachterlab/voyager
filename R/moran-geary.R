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
#' @param features Genes (\code{calculate*} SFE method and \code{run*}) or
#'   numeric columns of \code{colData(x)} (\code{colData*}) or any
#'   \code{\link{colGeometry}} (\code{colGeometryM*}) or
#'   \code{\link{annotGeometry}} (\code{annotGeometry*}) for which the
#'   univariate metric is to be computed. Default to \code{NULL}. When
#'   \code{NULL}, then the metric is computed for all genes with the values in
#'   the assay specified in the argument \code{exprs_values}. This can be
#'   parallelized with the argument \code{BPPARAM}.
#' @param exprs_values Integer scalar or string indicating which assay of x
#'   contains the expression values.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether
#'   and how computing the metric for numerous genes shall be parallelized.
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
#'   whose numeric columns of interest are to be used to compute the metric. Use
#'   \code{\link{colGeometryNames}} to look up names of the \code{sf} data
#'   frames associated with cells/spots.
#' @param annotGeometryName Name of a \code{annotGeometry} \code{sf} data frame
#'   whose numeric columns of interest are to be used to compute the metric. Use
#'   \code{\link{annotGeometryNames}} to look up names of the \code{sf} data
#'   frames associated with annotations.
#' @param sample_id Sample in the SFE object whose cells/spots to use.
#' @return For \code{calculate*}, a \code{DataFrame} with two columns: The
#'   first one is I for Moran's I or C for Geary's C, and the second one is K
#'   for sample kurtosis. For \code{run*}, a \code{SpatialFeatureExperiment}
#'   object with the Moran's I or Geary's C values added to a column of
#'   \code{rowData(x)}, whose name is specified in the \code{name} argument,
#'   with \code{sample_id} appended if applicable. For \code{colData},
#'   \code{colGeometry}, and \code{annotGeometry}, the results are added to the
#'   attributes of the corresponding columns, with regular names. For instance,
#'   Moran's I for nCounts for sample01 would in the "MoransI_sample01" of the
#'   nCounts column in \code{colData}.
#' @name calculateMoransI
#' @aliases calculateGearysC
#' @importFrom spdep moran geary Szero
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom SpatialFeatureExperiment SpatialFeatureExperiment
#' @importFrom SummarizedExperiment assay rowData<-
#' @importFrom SpatialFeatureExperiment colGraph annotGraph
#' @importFrom SingleCellExperiment colData rowData
NULL

.calc_univar_autocorr <- function(x, listw, fun, BPPARAM, returnDF = FALSE, ...) {
  if (is.null(listw))
    stop("The graph specified is absent from the SFE object.")
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  out_list <- bplapply(seq_len(nrow(x)), function(i) {
    fun(x[i,], listw, ...)
  }, BPPARAM = BPPARAM)
  if (returnDF) {
    out_list <- lapply(out_list, unlist, use.names = TRUE)
    out <- Reduce(rbind, out_list)
    if (!is.matrix(out)) out <- t(as.matrix(out))
    rownames(out) <- rownames(x)
    out <- DataFrame(out)
  } else {
    names(out_list) <- rownames(x)
    out <- out_list
  }
  return(out)
}

#' @rdname calculateMoransI
#' @export
setMethod("calculateMoransI", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL, returnDF = FALSE) {
  .calc_univar_autocorr(x, listw, fun = moran, BPPARAM = BPPARAM,
                        n = length(listw$neighbours), S0 = Szero(listw),
                        zero.policy = zero.policy, returnDF = returnDF)
})

#' @rdname calculateMoransI
#' @export
setMethod("calculateGearysC", "ANY", function(x, listw, BPPARAM = SerialParam(),
                                              zero.policy = NULL, returnDF = FALSE) {
  .calc_univar_autocorr(x, listw, fun = geary, BPPARAM = BPPARAM,
                        n = length(listw$neighbours),
                        n1 = length(listw$neighbours) - 1, S0 = Szero(listw),
                        zero.policy = zero.policy, returnDF = returnDF)
})
.check_sample_id <- SpatialFeatureExperiment:::.check_sample_id

.calc_univar_sfe_fun <- function(fun) {
  function(x, colGraphName, features = NULL, sample_id = NULL,
           exprs_values = "logcounts", BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    # Am I sure that I want to use logcounts as the default?
    sample_id <- .check_sample_id(x, sample_id)
    features <- .check_features(x, features)[["assay"]]
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
  mat <- t(as.matrix(df[, features, drop = FALSE]))
  if (anyNA(mat)) {
    stop("Only numeric columns without NA (within the sample_id) can be used.")
  }
  fun(mat, listw, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
}

.coldata_univar_fun <- function(fun) {
  name <- gsub("calculate", "", deparse(substitute(fun)))
  function(x, colGraphName, features, sample_id = NULL, BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    res <- .df_univar_autocorr(colData(x)[colData(x)$sample_id == sample_id,],
                               listw_use, features, fun,
                               BPPARAM, zero.policy, ...)
    name_attr <- paste(name, sample_id, sep = "_")
    for (i in seq_along(features)) {
      attr(colData(x)[[features[i]]], name_attr) <- res[[i]]
    }
    x
  }
}

#' @rdname calculateMoransI
#' @export
colDataMoransI <- .coldata_univar_fun(calculateMoransI)

#' @rdname calculateMoransI
#' @export
colDataGearysC <- .coldata_univar_fun(calculateGearysC)

.colgeom_univar_fun <- function(fun) {
  name <- gsub("calculate", "", deparse(substitute(fun)))
  function(x, colGeometryName, colGraphName, features, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    cg <- colGeometry(x, type = colGeometryName, sample_id = sample_id)
    res <- .df_univar_autocorr(cg, listw_use, features, fun, BPPARAM,
                               zero.policy, ...)
    name_attr <- paste(name, sample_id, sep = "_")
    cg_full <- colGeometry(x, type = colGeometryName)
    for (i in seq_along(features)) {
      attr(cg_full[[features[i]]], name_attr) <- res[[i]]
    }
    colGeometry(x, type = colGeometryName) <- cg_full
    x
  }
}

#' @rdname calculateMoransI
#' @export
colGeometryMoransI <- .colgeom_univar_fun(calculateMoransI)

#' @rdname calculateMoransI
#' @export
colGeometryGearysC <- .colgeom_univar_fun(calculateGearysC)

.annotgeom_univar_fun <- function(fun) {
  name <- gsub("calculate", "", deparse(substitute(fun)))
  function(x, annotGeometryName, annotGraphName, features, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- annotGraph(x, type = annotGraphName, sample_id = sample_id)
    .df_univar_autocorr(annotGeometry(x, type = annotGeometryName,
                                      sample_id = sample_id),
                        listw_use, features, fun, BPPARAM,
                        zero.policy, ...)
    name_attr <- paste(name, sample_id, sep = "_")
    ag_full <- annotGeometry(x, type = annotGeometryName)
    for (i in seq_along(features)) {
      attr(ag_full[[features[i]]], name_attr) <- res[[i]]
    }
    annotGeometry(x, type = annotGeometryName) <- ag_full
    x
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
             zero.policy, returnDF = TRUE)
  names(out)[1] <- name
  if (length(sampleIDs(x)) > 1L) {
    names(out) <- paste(names(out), sample_id, sep = "_")
  }
  rowData(x)[features, names(out)] <- out
  x
}

#' @rdname calculateMoransI
#' @export
runMoransI <- function(x, colGraphName, features, sample_id = NULL,
                       exprs_values = "logcounts", BPPARAM = SerialParam(),
                       zero.policy = NULL, name = "MoransI") {
  .sfe_univar_autocorr(x, colGraphName, features, sample_id, exprs_values,
                       calculateMoransI, BPPARAM, zero.policy, name = "MoransI")
}

#' @rdname calculateMoransI
#' @export
runGearysC <- function(x, colGraphName, features, sample_id = NULL,
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
                                              alternative = "greater",
                                              returnDF = FALSE, ...) {
  .calc_univar_autocorr(x, listw, fun = moran.mc, BPPARAM = BPPARAM,
                        nsim = nsim, zero.policy = zero.policy,
                        alternative = alternative, returnDF = returnDF, ...)
})

#' @rdname calculateMoranMC
#' @export
setMethod("calculateGearyMC", "ANY", function(x, listw, nsim,
                                              BPPARAM = SerialParam(),
                                              zero.policy = NULL,
                                              alternative = "greater",
                                              returnDF = FALSE, ...) {
  .calc_univar_autocorr(x, listw, fun = geary.mc, BPPARAM = BPPARAM,
                        nsim = nsim, zero.policy = zero.policy,
                        alternative = alternative, returnDF = returnDF, ...)
})

# For the `nsim` and `alternative` arguments
.calc_univar_sfe_fun_mc <- function(fun) {
  function(x, colGraphName, features, sample_id = NULL, nsim,
           exprs_values = "logcounts",
           BPPARAM = SerialParam(), zero.policy = NULL,
           alternative = "greater", ...) {
    .calc_univar_sfe_fun(fun)(x, colGraphName, features, sample_id,
                              exprs_values, BPPARAM, zero.policy, nsim = nsim,
                              alternative = alternative, returnDF = TRUE, ...)
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
  function(x, colGeometryName, colGraphName, features, sample_id = NULL, nsim,
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
  function(x, colGeometryName, colGraphName, features, sample_id = NULL, nsim,
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
  function(x, annotGeometryName, annotGraphName, features, sample_id = NULL, nsim,
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
             zero.policy, alternative, returnDF = TRUE, ...)
  # Convert results to DFrame. I'll write my own plotting function based on ggplot2.
  out <- lapply(out, function(o) {
    o$res <- I(list(o$res))
    DataFrame(unclass(o))
  })
  rns <- names(out)
  out <- Reduce(rbind, out)
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
runMoranMC <- function(x, colGraphName, features, sample_id = NULL, nsim,
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
runGearyMC <- function(x, colGraphName, features, sample_id = NULL, nsim,
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
          function(x, colGraphName, features, sample_id = NULL, order = 1,
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
                               sample_id = NULL, order = 1, method = "I",
                               BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .colgeom_univar_fun(calculateCorrelogram)(
    x, colGeometryName, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
colDataCorrelogram <- function(x, colGraphName, features, sample_id = NULL,
                               order = 1, method = "I", BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .coldata_univar_fun(calculateCorrelogram)(
    x, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
annotGeometryCorrelogram <- function(x, annotGeometryName, annotGraphName,
                                     features, sample_id = NULL, order = 1,
                                     method = "I", BPPARAM = SerialParam(),
                                     zero.policy = NULL, ...) {
  .annotgeom_univar_fun(calculateCorrelogram)(
    x, annotGeometryName, annotGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
runCorrelogram <- function(x, colGraphName, features, sample_id = NULL, order = 1,
                           method = "I", exprs_values = "logcounts",
                           BPPARAM = SerialParam(), zero.policy = NULL,
                           name = paste("Correlogram", method, sep = "_"), ...) {
  out <- calculateCorrelogram(x, colGraphName, features, sample_id, order,
                              method, exprs_values, BPPARAM, zero.policy, ...)
  if (method == "I") {
    out <- lapply(out, function(o) {
      colnames(o$res) <- c("I", "expectation", "variance")
      o
    })
  }
  out <- lapply(out, function(o) o$res)
  out_df <- DataFrame(res = I(out), row.names = features)
  if (length(sampleIDs(x)) > 1L) {
    name <- paste(name, sample_id, sep = "_")
  }
  names(out_df) <- name
  rowData(x)[features, names(out_df)] <- out_df
  x
}

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
#' @return For \code{calculateMoranPlot} and the \code{colData},
#'   \code{colGeometry}, and \code{annotGeometry} versions, a list of data
#'   frames that are the output of \code{\link{moran.plot}}. The names of the
#'   list are the names of the features. For \code{runMoranPlot}, the list is
#'   added to a column in \code{rowData(x)}. The plot is not made.
NULL

#' @rdname calculateMoranPlot
#' @export
setMethod("calculateMoranPlot", "ANY",
          function(x, listw, BPPARAM = SerialParam(),
                   zero.policy = NULL, ...) {
  .calc_univar_autocorr(x, listw, fun = moran.plot, BPPARAM, returnDF = FALSE,
                        plot = FALSE, return_df = TRUE, ...)
})

#' @rdname calculateMoranPlot
#' @export
setMethod("calculateMoranPlot", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateMoranPlot))

#' @rdname calculateMoranPlot
#' @export
colDataMoranPlot <- .coldata_univar_fun(calculateMoranPlot)

#' @rdname calculateMoranPlot
#' @export
colGeometryMoranPlot <- .colgeom_univar_fun(calculateMoranPlot)

#' @rdname calculateMoranPlot
#' @export
annotGeometryMoranPlot <- .annotgeom_univar_fun(calculateMoranPlot)

#' @rdname calculateMoranPlot
#' @export
runMoranPlot <- function(x, colGraphName, features, sample_id = NULL,
                         exprs_values = "logcounts", BPPARAM = SerialParam(),
                         zero.policy = NULL, name = "MoranPlot", ...) {
  out <- calculateMoranPlot(x, colGraphName, features, sample_id, exprs_values,
                            BPPARAM, zero.policy, ...)
  out_df <- DataFrame(res = I(out), row.names = features)
  if (length(sampleIDs(x)) > 1L) {
    name <- paste(name, sample_id, sep = "_")
  }
  names(out_df) <- name
  rowData(x)[features, names(out_df)] <- out_df
  x
}

#' Find clusters on the Moran plot
#'
#' The Moran plot plots the value at each location on the x axis, and the average
#' of the neighbors of each locations on the y axis. Sometimes clusters can be
#' seen on the Moran plot, indicating different types of neighborhoods.
#'
#' @inheritParams bluster::clusterRows
#' @inheritParams calculateMoranPlot
#' @param x A \code{SpatialFeatureExperiment} object with Moran plot computed
#' for the feature of interest. If the Moran plot for that feature has not been
#' computed for that feature in this sample_id, it will be calculated and stored
#' in \code{rowData}. See \code{\link{calculateMoranPlot}}.
#' @param feature One feature whose Moran plot is to be clustered.
#' @return The clusters are added to \code{colData} of the SFE object and the
#' SFE object is returned, just like in Seurat's \code{FindClusters}.
