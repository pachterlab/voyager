#' Bivariate spatial statistics
#'
#' These functions perform bivariate spatial analysis. In this version, the
#' bivariate global method supported are \code{\link{lee}},
#' \code{\link{lee.mc}}, and \code{\link{lee.test}} from \code{spdep}, and cross
#' variograms from \code{gstat} (use \code{cross_variogram} and
#' \code{cross_variogram_map} for \code{type} argument, see
#' \code{\link{variogram-internal}}). Global Lee statistic is computed by my own
#' implementation that is much faster than that in \code{spdep}. Bivariate local
#' methods supported are \code{\link{lee}} (use \code{locallee} for \code{type}
#' argument) and \code{\link{localmoran_bv}} a bivariate version of Local Moran
#' in \code{spdep}.
#'
#' @inheritParams calculateUnivariate
#' @param x A numeric matrix whose rows are features/genes, or a numeric vector
#'   (then \code{y} must be specified), or a \code{SpatialFeatureExperiment}
#'   (SFE) object with such a matrix in an assay.
#' @param y A numeric matrix whose rows are features/genes, or a numeric vector.
#'   Bivariate statics will be computed for all pairwise combinations of row
#'   names of x and row names of y, except in cross variogram where combinations
#'   within x and y are also computed.
#' @param feature1 ID or symbol of the first genes in SFE object, for the
#'   argument \code{x}.
#' @param feature2 ID or symbol of the second genes in SFE object, for the
#'   argument \code{x}. Mandatory if length of \code{feature1} is 1.
#' @return The \code{calculateBivariate} function returns a correlation matrix
#'   for global Lee, and the results for the each pair of genes for other
#'   methods. Global results are not stored in the SFE object. Some methods
#'   return one result for each pair of genes, while some return pairwise
#'   results for more than 2 genes jointly. Local results are stored in the
#'   \code{\link{localResults}} field in the SFE object, with name the
#'   concatenation the two gene names separated by two underscores (\code{__}).
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' library(scran)
#' library(SpatialFeatureExperiment)
#' library(SpatialExperiment)
#' sfe <- McKellarMuscleData()
#' sfe <- sfe[,sfe$in_tissue]
#' sfe <- logNormCounts(sfe)
#' gs <- modelGeneVar(sfe)
#' hvgs <- getTopHVGs(gs, fdr.threshold = 0.01)
#' g <- colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#'
#' # Matrix method
#' mat <- logcounts(sfe)[hvgs[1:5],]
#' df <- df2sf(spatialCoords(sfe), spatialCoordsNames(sfe))
#' out <- calculateBivariate(mat, type = "lee", listw = g)
#' out <- calculateBivariate(mat, type = "cross_variogram", coords_df = df)
#'
#' # SFE method
#' out <- calculateBivariate(sfe, type = "lee",
#' feature1 = c("Myh1", "Myh2", "Csrp3"), swap_rownames = "symbol")
#' out2 <- calculateBivariate(sfe, type = "lee.test", feature1 = "Myh1",
#' feature2 = "Myh2", swap_rownames = "symbol")
#' sfe <- runBivariate(sfe, type = "locallee", feature1 = "Myh1",
#' feature2 = "Myh2", swap_rownames = "symbol")
#' @name calculateBivariate
NULL

.to_mat1 <- function(x, name = "x") {
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    else if (ncol(x) == 1L) x <- t(x)
    if (is.null(rownames(x)) && nrow(x) == 1L) rownames(x) <- name
    x
}

.call_fun <- function(x, y, type, listw, zero.policy, coords_df,
                      simplify = TRUE, ...) {
    out <- if (use_graph(type)) fun(type)(x, y, listw, zero.policy = zero.policy, ...)
    else fun(type)(x, y, coords_df, ...)
    if (!simplify) out <- list(x__y = out)
    out
}

.call_fun_grid <- function(x, y, type, listw, zero.policy, coords_df,
                           BPPARAM, ...) {
    combs <- as.matrix(expand.grid(rownames(x), rownames(y)))
    out <- bplapply(seq_len(nrow(combs)), function(i) {
        .call_fun(x[combs[i,1],], y[combs[i,2],], type, listw, zero.policy,
                  coords_df, ...)
    }, BPPARAM = BPPARAM)
    names(out) <- paste(combs[,1], combs[,2], sep = "__")
    out
}

#' @rdname calculateBivariate
#' @export
setMethod("calculateBivariate", "ANY",
          function(x, y = NULL, type, listw = NULL, coords_df = NULL, BPPARAM = SerialParam(),
                   zero.policy = NULL, returnDF = TRUE, p.adjust.method = "BH",
                   name = NULL, ...) {
              if (is.character(type)) type <- get(type, mode = "S4")
              # x and y are expected to have genes in rows if they're matrices
              if (is.null(y)) {
                  if (is.vector(x) || min(dim(x)) == 1L) {
                      stop("y must be specified for vector x.")
                  } else {
                      x <- .to_mat1(x)
                      if (use_matrix(type)) {
                          out <- .call_fun(x, y = NULL, type = type,
                                           listw = listw,
                                           zero.policy = zero.policy,
                                           coords_df = coords_df, ...)
                      } else {
                          out <- .call_fun_grid(x, y = x, type, listw,
                                                zero.policy, coords_df, BPPARAM,
                                                ...)
                      }
                  }
              } else if (use_matrix(type) || (is.vector(x) && is.vector(y))) {
                  if (use_matrix(type)) {
                      x <- .to_mat1(x, "x")
                      y <- .to_mat1(y, "y")
                      is_diff_cells <- ncol(x) != ncol(y)
                  } else {
                      is_diff_cells <- length(x) != length(y)
                  }
                  if (is_diff_cells)
                      stop("x and y must have the same number of observations (cells).")
                  out <- .call_fun(x, y, type, listw = listw,
                                   zero.policy = zero.policy,
                                   coords_df = coords_df,
                                   simplify = !(is.vector(x) && is.vector(y)), ...)
              } else {
                  # fun only takes vectors for x and y and at least one of x and y is a matrix
                  x <- .to_mat1(x, "x")
                  y <- .to_mat1(y, "y")
                  if (ncol(x) != ncol(y)) {
                      stop("x and y must have the same number of observations (cells).")
                  }
                  if (is.null(rownames(x)) || is.null(rownames(y))) {
                      stop("Matrices x and y must have row names.")
                  }
                  out <- .call_fun_grid(x, y, type, listw, zero.policy,
                                        coords_df, BPPARAM, ...)
              }
              if (returnDF) {
                  if (is_local(type)) {
                      out <- reorganize_fun(type)(out, nb = listw$neighbours,
                                                  p.adjust.method = p.adjust.method)
                      out <- .value2df(out, use_geometry = FALSE)
                  } else {
                      if (is.null(name)) name <- info(type, "name")
                      out <- reorganize_fun(type)(out, name = name, ...)
                  }
              }
              if (length(out) == 1L && !is(out, "DataFrame")) out <- out[[1]]
              out
          })

#' @rdname calculateBivariate
#' @export
setMethod("calculateBivariate", "SpatialFeatureExperiment",
          # For matrix, specifically for global Lee
          function(x, type, feature1, feature2 = NULL, colGraphName = 1L,
                   colGeometryName = 1L, sample_id = "all",
                   exprs_values = "logcounts", BPPARAM = SerialParam(),
                   zero.policy = NULL, returnDF = TRUE, p.adjust.method = "BH",
                   swap_rownames = NULL, name = NULL, ...) {
              sample_id <- .check_sample_id(x, sample_id, one = FALSE)
              if (is.character(type)) type <- get(type, mode = "S4")
              if (is.null(feature2) && length(feature1) == 1L) {
                  stop("feature2 must be specified when feature1 has length 1.")
              }
              out <- lapply(sample_id, function(s) {
                  feature1 <- .check_features(x, feature1, swap_rownames = swap_rownames)[["assay"]]
                  mat1 <- assay(x, exprs_values)[feature1, colData(x)$sample_id == s, drop = FALSE]
                  if (!is.null(feature2)) {
                      feature2 <- .check_features(x, feature2, swap_rownames = swap_rownames)[["assay"]]
                      mat2 <- assay(x, exprs_values)[feature2, colData(x)$sample_id == s, drop = FALSE]
                  } else mat2 <- NULL

                  if (use_graph(type)) {
                      listw_use <- colGraph(x, type = colGraphName, sample_id = s)
                      o <- calculateBivariate(mat1, mat2, listw = listw_use,
                                              type = type,
                                              BPPARAM = BPPARAM,
                                              zero.policy = zero.policy,
                                              returnDF = returnDF, p.adjust.method = p.adjust.method,
                                              name = name, ...
                      )
                  } else {
                      cg <- colGeometry(x, colGeometryName, sample_id = s)
                      cg <- .get_coords_df(x, cg, s, exprs_values, swap_rownames, ...)
                      o <- calculateBivariate(mat1, mat2, coords_df = cg,
                                              type = type, BPPARAM = BPPARAM,
                                              returnDF = returnDF,
                                              p.adjust.method = p.adjust.method,
                                              name = name, ...)
                  }
                  o
              })
              names(out) <- sample_id
              if (length(sample_id) == 1L) out <- out[[1]]
              out
          })

#' @rdname calculateBivariate
#' @export
runBivariate <- function(x, type, feature1, feature2 = NULL, colGraphName = 1L,
                         colGeometryName = 1L, sample_id = "all",
                         exprs_values = "logcounts", BPPARAM = SerialParam(),
                         swap_rownames = NULL,
                         zero.policy = NULL,
                         p.adjust.method = "BH", name = NULL, ...) {
    if (is.character(type)) type <- get(type, mode = "S4")
    if (!is_local(type)) {
        stop("Global bivariate results can't be stored in the SFE object.",
             " Use calculateBivariate instead.")
    }
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (is.null(name)) name <- info(type, "name")
    other_args <- list(...)
    if (use_graph(type))
        g <- colGraph(x, type = colGraphName, sample_id = sample_id[1])
    else g <- NULL
    params <- c(info(type, c("name", "package")),
                list(version = packageVersion(info(type, "package")),
                     zero.policy = zero.policy,
                     p.adjust.method = p.adjust.method,
                     graph_params = attr(g, "method")), other_args)

    old_params <- getParams(x, name, local = TRUE)
    .check_old_params(params, old_params, name, args_not_check(type))

    feature1_id <- .symbol2id(x, feature1, swap_rownames)
    if (!is.null(feature2)) feature2_id <- .symbol2id(x, feature2, swap_rownames)
    else feature2_id <- NULL
    for (s in sample_id) {
        out <- calculateBivariate(x, type, feature1_id, feature2_id,
                                  colGraphName, colGeometryName, s,
                                  exprs_values, BPPARAM, zero.policy,
                                  returnDF = TRUE,
                                  p.adjust.method = p.adjust.method, ...
        )
        if (!is.null(swap_rownames)) {
            if (is.null(feature2)) feature2 <- feature1
            comb <- expand.grid(feature1, feature2)
            feature_use <- paste(comb[,1], comb[,2], sep = "__")
            names(out) <- feature_use
        }
        x <- .add_localResults_info(x, sample_id = s,
                                    name = name, features = feature_use,
                                    res = out, params = params)
    }
    x
}
