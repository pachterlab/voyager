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
#' @return For \code{calculate*}, a \code{DataFrame} with two columns: The first
#'   one is I for Moran's I or C for Geary's C, and the second one is K for
#'   sample kurtosis. For \code{run*}, a \code{SpatialFeatureExperiment} object
#'   with the Moran's I or Geary's C values added to a column of
#'   \code{rowData(x)}, whose name is specified in the \code{name} argument,
#'   with \code{sample_id} appended if applicable. For \code{colData},
#'   \code{colGeometry}, and \code{annotGeometry}, the results are added to an
#'   attribute of the data frame called \code{featureData}, which is a DataFrame
#'   analogous to \code{rowData} for the gene count matrix. New column names in
#'   \code{featureData} would follow the same rules as in \code{rowData}. (I
#'   need to write many examples to make it clear to users.)
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
  if (is(df, "sf")) df <- st_drop_geometry(df)
  mat <- t(as.matrix(df[, features, drop = FALSE]))
  if (anyNA(mat)) {
    stop("Only numeric columns without NA (within the sample_id) can be used.")
  }
  fun(mat, listw, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
}

.MoransI2df <- function(out, name) {
  # out should already be a DFrame when this function is called
  names(out)[1] <- name
  out
}

.coldata_univar_fun <- function(fun, to_df_fun, name, to_df_params = list()) {
  function(x, colGraphName, features, sample_id = NULL, BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    res <- .df_univar_autocorr(colData(x)[colData(x)$sample_id == sample_id,],
                               listw_use, features, fun,
                               BPPARAM, zero.policy, ...)
    colData(x) <- .add_fd(x, colData(x), res, features, sample_id, to_df_fun,
                          name, to_df_params)
    x
  }
}

#' @rdname calculateMoransI
#' @export
colDataMoransI <- .coldata_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
colDataGearysC <- .coldata_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

.colgeom_univar_fun <- function(fun, to_df_fun, name, to_df_params = list()) {
  function(x, colGeometryName, colGraphName, features, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- colGraph(x, type = colGraphName, sample_id = sample_id)
    cg <- colGeometry(x, type = colGeometryName, sample_id = sample_id)
    res <- .df_univar_autocorr(cg, listw_use, features, fun, BPPARAM,
                               zero.policy, ...)
    colGeometry(x, colGeometryName) <- .add_fd(x, colGeometry(x, colGeometryName),
                                               res, features, sample_id,
                                               to_df_fun, name, to_df_params)
    x
  }
}

#' @rdname calculateMoransI
#' @export
colGeometryMoransI <- .colgeom_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
colGeometryGearysC <- .colgeom_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

.annotgeom_univar_fun <- function(fun, to_df_fun, name, to_df_params = list()) {
  function(x, annotGeometryName, annotGraphName, features, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id)
    listw_use <- annotGraph(x, type = annotGraphName, sample_id = sample_id)
    res <- .df_univar_autocorr(annotGeometry(x, type = annotGeometryName,
                                             sample_id = sample_id),
                               listw_use, features, fun, BPPARAM,
                               zero.policy, ...)
    annotGeometry(x, annotGeometryName) <- .add_fd(x, annotGeometry(x, annotGeometryName),
                                                   res, features, sample_id,
                                                   to_df_fun, name, to_df_params)
    x
  }
}

#' @rdname calculateMoransI
#' @export
annotGeometryMoransI <- .annotgeom_univar_fun(calculateMoransI, .MoransI2df, "MoransI")

#' @rdname calculateMoransI
#' @export
annotGeometryGearysC <- .annotgeom_univar_fun(calculateGearysC, .MoransI2df, "GearysC")

.sfe_univar_autocorr <- function(x, colGraphName, features, sample_id,
                                 exprs_values, fun, BPPARAM, zero.policy,
                                 name) {
  sample_id <- .check_sample_id(x, sample_id)
  out <- fun(x, colGraphName, features, sample_id, exprs_values, BPPARAM,
             zero.policy)
  out <- .MoransI2df(out, name)
  out <- .add_name_sample_id(x, out, sample_id)
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
#' is returned. For the colData, colGeometry, and annotGeometry versions, the
#' results are added to the \code{featureData} attribute of the data frame of
#' interest in a manner analogous to \code{rowData}.
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
  function(x, colGraphName, features, sample_id = NULL, nsim,
           exprs_values = "logcounts",
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

.coldata_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, colGraphName, features, sample_id = NULL, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .coldata_univar_fun(fun, to_df_fun, name)(x, colGraphName,
                                              features, sample_id, BPPARAM,
                                              zero.policy, nsim = nsim,
                                              alternative = alternative, ...)
  }
}

.MoranMC2df <- function(out, name) {
  # Convert results to DFrame. I'll write my own plotting function based on ggplot2.
  out <- lapply(out, function(o) {
    o$res <- I(list(o$res))
    DataFrame(unclass(o))
  })
  rns <- names(out)
  out <- Reduce(rbind, out)
  rownames(out) <- rns
  colnames(out) <- paste(name, colnames(out), sep = "_")
  out
}

#' @rdname calculateMoranMC
#' @export
colDataMoranMC <- .coldata_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
colDataGearyMC <- .coldata_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

.colgeom_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, colGeometryName, colGraphName, features, sample_id = NULL, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .colgeom_univar_fun(fun, to_df_fun, name)(x, colGeometryName, colGraphName,
                                              features, sample_id, BPPARAM,
                                              zero.policy, nsim = nsim,
                                              alternative = alternative, ...)
  }
}
#' @rdname calculateMoranMC
#' @export
colGeometryMoranMC <- .colgeom_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
colGeometryGearyMC <- .colgeom_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

.annotgeom_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, annotGeometryName, annotGraphName, features, sample_id = NULL, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .annotgeom_univar_fun(fun, to_df_fun, name)(
      x, annotGeometryName, annotGraphName, features, sample_id, BPPARAM,
      zero.policy, nsim = nsim, alternative = alternative, ...)
  }
}

#' @rdname calculateMoranMC
#' @export
annotGeometryMoranMC <- .annotgeom_univar_fun_mc(calculateMoranMC, .MoranMC2df, "MoranMC")

#' @rdname calculateMoranMC
#' @export
annotGeometryGearyMC <- .annotgeom_univar_fun_mc(calculateGearyMC, .MoranMC2df, "GearyMC")

.sfe_univar_mc <- function(x, colGraphName, features, sample_id, nsim,
                           exprs_values, fun, BPPARAM, zero.policy, alternative,
                           name, ...) {
  out <- fun(x, colGraphName, features, sample_id, nsim, exprs_values, BPPARAM,
             zero.policy, alternative, ...)
  out <- .MoranMC2df(out, name)
  out <- .add_name_sample_id(x, out, sample_id)
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
#' @return For \code{calculateCorrelogram}, a list of
#'   \code{spcor} objects, each element of which correslonds to a feature. For
#'   \code{runCorrelogram}, the \code{res} field of the \code{spcor} is taken
#'   and put in a list column in \code{rowData(x)}, and the SFE object with the
#'   new \code{rowData} is returned. For the colData, colGeometry, and
#'   annotGeometry versions, the results are added to an attribute of the data
#'   frame of interest called \code{featureData}, in a manner analogous to
#'   \code{rowData}.
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

.correlogram2df <- function(out, name, method) {
  if (method == "I") {
    out <- lapply(out, function(o) {
      colnames(o$res) <- c("I", "expectation", "variance")
      o
    })
  }
  out <- lapply(out, function(o) o$res)
  out_df <- DataFrame(res = I(out))
  names(out_df) <- name
  out_df
}

#' @rdname calculateCorrelogram
#' @export
colGeometryCorrelogram <- function(x, colGeometryName, colGraphName, features,
                               sample_id = NULL, order = 1, method = "I",
                               BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .colgeom_univar_fun(calculateCorrelogram, .correlogram2df,
                      name = paste("Correlogram", method, sep = "_"),
                      to_df_params = list(method = method))(
    x, colGeometryName, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
colDataCorrelogram <- function(x, colGraphName, features, sample_id = NULL,
                               order = 1, method = "I", BPPARAM = SerialParam(),
                               zero.policy = NULL, ...) {
  .coldata_univar_fun(calculateCorrelogram, .correlogram2df,
                      name = paste("Correlogram", method, sep = "_"),
                      to_df_params = list(method = method))(
    x, colGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
annotGeometryCorrelogram <- function(x, annotGeometryName, annotGraphName,
                                     features, sample_id = NULL, order = 1,
                                     method = "I", BPPARAM = SerialParam(),
                                     zero.policy = NULL, ...) {
  .annotgeom_univar_fun(calculateCorrelogram, .correlogram2df,
                        name = paste("Correlogram", method, sep = "_"),
                        to_df_params = list(method = method))(
    x, annotGeometryName, annotGraphName, features, sample_id, BPPARAM = BPPARAM,
    zero.policy = zero.policy, order = order, method = method, ...)
}

#' @rdname calculateCorrelogram
#' @export
runCorrelogram <- function(x, colGraphName, features, sample_id = NULL, order = 1,
                           method = "I", exprs_values = "logcounts",
                           BPPARAM = SerialParam(), zero.policy = NULL,
                           name = paste("Correlogram", method, sep = "_"), ...) {
  sample_id <- .check_sample_id(x, sample_id)
  out <- calculateCorrelogram(x, colGraphName, features, sample_id, order,
                              method, exprs_values, BPPARAM, zero.policy, ...)
  out <- .correlogram2df(out, name, method)
  rownames(out) <- features
  out <- .add_name_sample_id(x, out, sample_id)
  rowData(x)[features, names(out)] <- out
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
#' @return For \code{calculateMoranPlot}, a list of data frames that are the
#'   output of \code{\link{moran.plot}}. The names of the list are the names of
#'   the features. For \code{runMoranPlot}, the list is added to a column in
#'   \code{rowData(x)}. For the colData, colGeometry, and annotGeometry
#'   versions, the results are added to an attribute of the data frame of
#'   interest called \code{featureData}, in a manner analogous to
#'   \code{rowData}. The plot is not made.
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
runMoranPlot <- function(x, colGraphName, features, sample_id = NULL,
                         exprs_values = "logcounts", BPPARAM = SerialParam(),
                         zero.policy = NULL, name = "MoranPlot", ...) {
  sample_id <- .check_sample_id(x, sample_id)
  out <- calculateMoranPlot(x, colGraphName, features, sample_id, exprs_values,
                            BPPARAM, zero.policy, ...)
  out <- .MoranPlot2df(out, name)
  rownames(out) <- features
  out <- .add_name_sample_id(x, out, sample_id)
  rowData(x)[features, names(out)] <- out
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
#' @param colGeometryName Name of colGeometry from which the feature is from. If
#' specified, then this function will only look in that colGeometry.
#' @param annotGeometryName Name of annotGeometry from which the feature is from.
#' If both \code{colGeometryName} and \code{annotGeometryName} are specified,
#' then a warning is issued and \code{colGeometryName} will be used. If neither
#' is specified, then this function will only look in \code{rownames(x)} and
#' \code{colData(x)}.
#' @param features Features whose Moran plot are to be cluster. Features whose
#' Moran plots have not been computed will be skipped, with a warning.
#' @return A data frame each column of which is a factor for cluster membership
#' of each feature. The column names are the features.
#' @importFrom bluster clusterRows
#' @export
clusterMoranPlot <- function(x, features, BLUSPARAM, sample_id = NULL,
                             name = "MoranPlot",
                             colGeometryName = NULL,
                             annotGeometryName = NULL) {
  sample_id <- .check_sample_id(x, sample_id)
  colname_use <- paste(name, sample_id, sep = "_")
  if (is.null(colGeometryName) && is.null(annotGeometryName)) {
    features <- .check_features(x, features)
    if (length(features[["assay"]]) && colname_use %in% names(rowData(x))) {
      mps_assay <- rowData(x)[features[["assay"]], colname_use]
    } else mps_assay <- NULL
    fd <- attr(colData(x), "featureData")
    if (is.null(fd)) mps_coldata <- NULL
    else if (length(features[["coldata"]]) && colname_use %in% names(fd)) {
      mps_coldata <- fd[features[["coldata"]], colname_use]
    } else mps_coldata <- NULL
    mps <- c(mps_assay, mps_coldata)
    if (!is.null(mps))
      names(mps) <- c(features[["assay"]], features[["coldata"]])
  } else if (!is.null(colGeometryName)) {
    if (!is.null(annotGeometryName)) {
      warning("Using colGeometryName rather than annotGeometryName.")
    }
    fd <- attr(colGeometry(x, colGeometryName, sample_id), "featureData")
    mps <- fd[features, colname_use]
    names(mps) <- features
  } else if (!is.null(annotGeometryName)) {
    fd <- attr(annotGeometry(x, annotGeometryName, sample_id), "featureData")
    mps <- fd[features, colname_use]
    names(mps) <- features
  }
  # isTRUE because is.na can return value of length > 1
  na_inds <- vapply(mps, function(t) isTRUE(is.na(t)), FUN.VALUE = logical(1))
  if (all(na_inds))
    stop("None of the features requested have Moran plots computed.")
  if (any(na_inds))
    warning("Skipping features that don't have Moran plots computed.")
  mps <- mps[!na_inds]
  out <- lapply(mps, function(mp) clusterRows(mp[,c("x", "wx")], BLUSPARAM))
  as.data.frame(out, row.names = colnames(x)[colData(x)$sample_id == sample_id])
}
