# Multivariate
# 1. GWPCA, from GWmodels
# 3. I really wish that MSPA can run faster. If I really want to use it, I'll reimplement it
# 4. Moran eigenmaps and spatial filtering
# 5. Spatially informed clustering


#' Multivariate spatial data analysis
#'
#' These functions perform multivariate spatial data analysis, usually spatially
#' informed dimension reduction.
#'
#' @inheritParams calculateUnivariate
#' @inheritParams scater::runPCA
#' @param sample_action Character, either "joint" or "separate". Spatial methods
#'   depend on the spatial coordinates and/or spatial neighborhood graph, which
#'   is why \code{SpatialExperiment} uses \code{sample_id} to keep coordinates
#'   from different samples separate. Some spatial methods can be sensibly run
#'   jointly for multiple samples. In this case, "joint" will run the method
#'   jointly for all samples, and "separate" will run the method separately for
#'   each sample and concatenate the results. If the output of each sample is a
#'   3D array, as in GWPCA, then the \code{abind} package needs to be installed
#'   to concatenate the results of the different samples.
#' @param dest Character, either "reducedDim" or "colData". If the output of the
#'   multivariate method is a matrix or array, as in spatially informed
#'   dimension reduction, then the only option is "reducedDim", so the results
#'   will be stored in \code{\link{reducedDim}} of the SFE object. If the output
#'   is a vector, as in the multivariate version of \code{\link{localC}}, then
#'   choosing "colData" will store the output in a new column in
#'   \code{colData(x)}.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether
#'   and how computing the metric for numerous genes shall be parallelized. This
#'   is to parallelize computation across multiple samples when there are a
#'   large number of samples. Be cautious if using an optimized BLAS for matrix
#'   operations that supports multithreading.
#' @param ... Extra arguments passed to the specific multivariate method. For
#'   example, see \code{\link{multispati_rsp}} for arguments for MULTISPATI PCA.
#' @return In \code{calculateMultivariate}, a matrix for cell embeddings whose
#'   attributes include loadings and eigenvalues if relevant, ready to be added
#'   to the SFE object with \code{reducedDim} setter. For \code{run*}, a
#'   \code{SpatialFeatureExperiment} object with the results added. See Details
#'   for where the results are stored.
#' @export
#' @name calculateMultivariate
#' @examples
#' # example code
#' library(SFEData)
#' library(scater)
#' library(scran)
#' sfe <- McKellarMuscleData()
#' sfe <- logNormCounts(sfe)
#' gvs <- modelGeneVar(sfe)
#' hvgs <- getTopHVGs(gvs, fdr.threshold = 0.05)
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' sfe <- runMultivariate(sfe, "multispati", subset_row = hvgs)
NULL

#' @rdname calculateMultivariate
#' @export
setMethod("calculateMultivariate", c("ANY", "SFEMethod"),
          function(x, type, listw = NULL, ...) {
              if (use_graph(type))
                  res <- fun(type)(x, listw, ...)
              else res <- fun(type)(x, ...)
              res <- reorganize_fun(type)(res)
              res
          })

#' @rdname calculateMultivariate
#' @export
setMethod("calculateMultivariate", c("ANY", "character"),
          function(x, type, listw = NULL, ...) {
              type <- get(type, mode = "S4")
              calculateMultivariate(x, type, listw, ...)
          })

# joint: MULTISPATI: Moran's I across all the samples, if spatial clustering
# across samples is needed. More relevant for biological replica and sequential sections?
# separate: Maximize Moran's I times variance explained in each sample separately.
# What if some genes have different spatial autocorrelation characteristics between
# case and control?

#' @rdname calculateMultivariate
#' @export
setMethod("calculateMultivariate", "SpatialFeatureExperiment",
          function(x, type, colGraphName = 1L,
                   ntop = 500, subset_row = NULL, exprs_values = "logcounts",
                   sample_action = c("joint", "separate"),
                   BPPARAM = SerialParam(), ...) {
              sample_id <- sampleIDs(x)
              if (is.character(type)) type <- get(type, mode = "S4")
              sample_action <- match.arg(sample_action)
              if (!is_joint(type) && sample_action == "joint")
                  sample_action <- "separate"

              mat <- assay(sfe, exprs_values)
              o <- scater:::.get_mat_for_reddim(mat, subset_row=subset_row,
                                                ntop=ntop,
                                                scale=FALSE, get.var=TRUE)
              mat <- o$x
              cv <- o$v
              if (length(sample_id) > 1L) {
                  bcs <- colnames(x)
                  if (use_graph(type))
                      listws <- colGraphs(x, sample_id = sample_id,
                                          name = colGraphName)
                  # outputs must have rownames or names
                  # I'll check if localC output has names
                  if (sample_action == "joint") {
                      if (use_graph(type)) {
                          W <- multi_listw2sparse(listws)
                          W <- W[bcs, bcs]
                      } else W <- NULL
                      out <- calculateMultivariate(mat, type, W, ...)
                  } else {
                      out <- bplapply(seq_along(sample_id), function(i) {
                          m <- mat[, x$sample_id == sample_id[i]]
                          l <- if (use_graph(type)) listws[[i]] else NULL
                          calculateMultivariate(m, type, l, ...)
                      }, BPPARAM = BPPARAM)

                      # Deal with attributes like loadings and eigenvalues from
                      # multiple samples
                      extra_attrs <- setdiff(names(attributes(out[[1]])),
                                             c("dim", "dimnames"))
                      # Need to double check gwpca output
                      attrs_combined <- list()
                      if (length(extra_attrs)) {
                          # Stack matrices into 3D arrays, vectors into matrices
                          for (e in extra_attrs) {
                              es <- lapply(out, function(o) attr(o, e))
                              if (is.vector(es[[1]]) && is.atomic(es[[1]]) &&
                                  length(unique(lengths(es))) == 1L) {
                                  orig_names <- names(es[[1]])
                                  comb <- do.call(cbind, es)
                                  colnames(comb) <- sample_id
                              } else if (is.matrix(es[[1]]) &&
                                         length(unique(lapply(es, dim))) == 1L) {
                                  comb <- simplify2array(es)
                                  dimnames(comb) <- c(dimnames(es[[1]]), list(sample_id))
                              } else {
                                  comb <- setNames(es, sample_id)
                              }
                              attrs_combined <- c(attrs_combined, comb)
                          }
                          names(attrs_combined) <- extra_attrs
                      }

                      is_vector <- all(vapply(out, is.vector, FUN.VALUE = logical(1)))
                      is_matrix <- all(vapply(out, is.matrix, FUN.VALUE = logical(1)))
                      is_array <- all(vapply(out, function(o) is.array(o) & length(dim(o)) > 2L,
                                             FUN.VALUE = logical(1)))
                      if (is_vector) {
                          out <- unlist(out)
                          out <- out[bcs]
                      } else if (is_matrix) {
                          out <- do.call(rbind, out)
                          out <- out[bcs,]
                      } else if (is_array) {
                          # For gwpca, need to check rownames
                          rlang::check_installed("abind")
                          out <- do.call(abind::abind, out)
                          out <- out[bcs,,]
                      }
                      attributes(out) <- c(attributes(out), attrs_combined)
                  }
              } else {
                  listw <- if (use_graph(type)) colGraph(x, colGraphName, sample_id) else NULL
                  out <- calculateMultivariate(mat, type, listw, ...)
              }
              out
          })

#' @rdname calculateMultivariate
#' @export
runMultivariate <- function(x, type, colGraphName = 1L,
                            ntop = 500, subset_row = NULL, exprs_values = "logcounts",
                            sample_action = c("joint", "separate"),
                            BPPARAM = SerialParam(), name = NULL,
                            dest = c("reducedDim", "colData"), ...) {
    dest <- match.arg(dest)
    if (is.character(type)) type <- get(type, mode = "S4")
    if (is.null(name)) name <- info(type, "name")
    out <- calculateMultivariate(x, type, colGraphName = colGraphName,
                                 ntop = ntop, subset_row = subset_row, exprs_values = exprs_values,
                                 sample_action = sample_action,
                                 BPPARAM = BPPARAM, name = name, ...)
    if (is.array(out) && dest == "colData") {
        message("Matrix or array outputs can only be stored in reducedDims.")
        dest <- "reducedDim"
    }
    if (dest == "reducedDim") {
        reducedDim(x, name) <- out
    } else {
        colData(x)[[name]] <- out
    }
    x
}
