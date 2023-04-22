# Multivariate
# 3. I really wish that MSPA can run faster. If I really want to use it, I'll reimplement it
# 4. Moran eigenmaps and spatial filtering
# 5. Spatially informed clustering
# GWPCA: I haven't completely given up. It's not that interesting in the Visium
# muscle data. It might be interesting in some smFISH dataset. Then I might need
# to reimplement it myself using Irlba or RSpectra and parallel programming to
# make it scalable enough to be more useful. I don't have enough time for that
# for Bioc 3.17. I would also need to write my own geom for the loading glyph
# plot to only show genes with top loadings rather than all genes.
# For both MSPA and GWPCA: the interesting part is not the cell embeddings, but
# the results about genes. I may write a different function or in whichever way
# note that the results will go into rowData.

#' Multivariate spatial data analysis
#'
#' These functions perform multivariate spatial data analysis, usually spatially
#' informed dimension reduction.
#'
#' For the argument \code{type}, this package supports "multispati" for
#' MULTISPATI PCA, "localC_multi" for a multivariate generalization of Geary's
#' C, "localC_perm_multi" for the multivariate Geary's C with permutation
#' testing, and "gwpca" for geographically weighted PCA.
#'
#' @inheritParams calculateUnivariate
#' @inheritParams scater::runPCA
#' @param sample_action Character, either "joint" or "separate". Spatial methods
#'   depend on the spatial coordinates and/or spatial neighborhood graph, which
#'   is why \code{SpatialExperiment} uses \code{sample_id} to keep coordinates
#'   from different samples separate. Some spatial methods can be sensibly run
#'   jointly for multiple samples. In this case, "joint" will run the method
#'   jointly for all samples, and "separate" will run the method separately for
#'   each sample and concatenate the results.
#' @param dest Character, either "reducedDim" or "colData". If the output of the
#'   multivariate method is a matrix or array, as in spatially informed
#'   dimension reduction, then the only option is "reducedDim", so the results
#'   will be stored in \code{\link{reducedDim}} of the SFE object. If the output
#'   is a vector, as in the multivariate version of \code{\link{localC}}, then
#'   it will be sotred in \code{colData}. Data frame output, such as from
#'   \code{localC_perm}, can be stored in either \code{reducedDim} or
#'   \code{colData}.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying whether
#'   and how computing the metric for numerous genes shall be parallelized. This
#'   is to parallelize computation across multiple samples when there are a
#'   large number of samples. Be cautious if using an optimized BLAS for matrix
#'   operations that supports multithreading.
#' @param transposed Logical, whether the matrix has genes in columns and cells
#'   in rows.
#' @param ... Extra arguments passed to the specific multivariate method. For
#'   example, see \code{\link{multispati_rsp}} for arguments for MULTISPATI PCA.
#'   See \code{\link{localC}} for arguments for "localC_multi" and
#'   "localC_perm_multi".
#' @return In \code{calculateMultivariate}, a matrix for cell embeddings whose
#'   attributes include loadings and eigenvalues if relevant, ready to be added
#'   to the SFE object with \code{reducedDim} setter. For \code{run*}, a
#'   \code{SpatialFeatureExperiment} object with the results added. See Details
#'   for where the results are stored.
#' @references Dray, S., Said, S. and Debias, F. (2008) Spatial ordination of
#' vegetation data using a generalization of Wartenberg's multivariate spatial
#' correlation. Journal of vegetation science, 19, 45-56.
#'
#' Anselin, L. (2019), A Local Indicator of Multivariate Spatial Association:
#' Extending Geary's c. Geogr Anal, 51: 133-150. doi:10.1111/gean.12164
#'
#' @export
#' @importFrom SpatialFeatureExperiment colGraphs
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
          function(x, type, listw = NULL, transposed = FALSE, zero.policy = TRUE,
                   p.adjust.method = "BH", ...) {
              if (info(type, "variate") != "multi")
                  stop("`type` must be a multivariate method.")
              if (!info(type, "package") %in% c("spdep", "Voyager", NA)) {
                  rlang::check_installed(info(type, "package"))
              }
              # x is expect to have genes in rows and cells in columns
              if (!transposed) x <- t(x)
              if (use_graph(type))
                  res <- fun(type)(x, listw, ...)
              else res <- fun(type)(x, ...)
              if (is_local(type)) {
                  res <- reorganize_fun(type)(res, nb = listw$neighbours,
                                              p.adjust.method = p.adjust.method)
              } else res <- reorganize_fun(type)(res)
              res
          })

#' @rdname calculateMultivariate
#' @export
setMethod("calculateMultivariate", c("ANY", "character"),
          function(x, type, listw = NULL, transposed = FALSE, ...) {
              type <- get(type, mode = "S4")
              calculateMultivariate(x, type, listw, transposed, ...)
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
                   subset_row = NULL, exprs_values = "logcounts",
                   sample_action = c("joint", "separate"),
                   BPPARAM = SerialParam(), ...) {
              sample_id <- sampleIDs(x)
              if (is.character(type)) type <- get(type, mode = "S4")
              sample_action <- match.arg(sample_action)
              if (!is_joint(type) && sample_action == "joint")
                  sample_action <- "separate"

              mat <- assay(x, exprs_values)
              if (!is.null(subset_row)) mat <- mat[subset_row,]
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
                                  dimnames(comb)[[3]] <- sample_id
                              } else {
                                  comb <- setNames(es, sample_id)
                              }
                              attrs_combined <- c(attrs_combined, list(comb))
                          }
                          names(attrs_combined) <- extra_attrs
                      }

                      is_vector <- all(vapply(out, is.vector, FUN.VALUE = logical(1)))
                      is_matrix <- all(vapply(out, is.matrix, FUN.VALUE = logical(1)))
                      is_df <- all(vapply(out, function(o) is.data.frame(o) | is(o, "DataFrame"),
                                          FUN.VALUE = logical(1)))
                      if (is_vector) {
                          out <- unlist(out)
                          out <- out[bcs]
                      } else if (is_matrix || is_df) {
                          out <- do.call(rbind, out)
                          out <- out[bcs,]
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
                            subset_row = NULL, exprs_values = "logcounts",
                            sample_action = c("joint", "separate"),
                            BPPARAM = SerialParam(), name = NULL,
                            dest = c("reducedDim", "colData"), ...) {
    dest <- match.arg(dest)
    if (is.character(type)) type <- get(type, mode = "S4")
    if (is.null(name)) name <- info(type, "name")
    out <- calculateMultivariate(x, type, colGraphName = colGraphName,
                                 subset_row = subset_row, exprs_values = exprs_values,
                                 sample_action = sample_action,
                                 BPPARAM = BPPARAM, ...)
    if (is.array(out) && dest == "colData") {
        message("Matrix or array outputs can only be stored in reducedDims.")
        dest <- "reducedDim"
    }
    if (is.vector(out) && dest == "reducedDim") {
        message("Vector output can only be stored in colData.")
        dest <- "colData"
    }
    if (dest == "reducedDim") {
        rownames(out) <- colnames(x)
        reducedDim(x, name) <- out
    } else {
        if (is.vector(out)) colData(x)[[name]] <- out
        if (is.data.frame(out) || is(out, "DataFrame")) {
            inds <- !grepl(name, names(out))
            names(out)[inds] <- paste(name, names(out)[inds], sep = "_")
            colData(x) <- cbind(colData(x), out)
        }
    }
    x
}
