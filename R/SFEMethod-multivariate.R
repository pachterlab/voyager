# SFEMethod objects for multivariate methods

# MULTISPATI PCA==================

#' A faster implementation of MULTISPATI PCA
#'
#' This implementation uses the \code{RSpectra} package to efficiently compute a
#' small subset of eigenvalues and eigenvectors, as a small subset is typically
#' used. Hence it's much faster and memory efficient than the original
#' implementation in \code{adespatial}. However, this implementation here does
#' not support row and column weighting other than the standard ones for PCA.,
#' so the \code{adespatial} implementation is more general.
#'
#' @param x A matrix whose columns are features and rows are cells.
#' @param listw A \code{listw} object, a spatial neighborhood graph for the
#'   cells in \code{x}. The length must be equal to the number of row of
#'   \code{x}.
#' @param nfposi Number of positive eigenvalues and their eigenvectors to
#'   compute.
#' @param nfnega Number of nega eigenvalues and their eigenvectors to compute.
#'   These indicate negative spatial autocorrelation.
#' @param scale Logical, whether to scale the data.
#' @return A matrix for the cell embeddings in each spatial PC, with attribute
#'   \code{loading} for the eigenvectors or gene loadings, and attribute
#'   \code{eig} for the eigenvalues.
#' @note
#' Eigen decomposition will fail if any feature has variance zero leading to NaN
#' in the scaled matrix.
#' @references
#' Dray, S., Said, S. and Debias, F. (2008) Spatial ordination of vegetation data using a generalization of Wartenberg's multivariate spatial correlation. Journal of vegetation science, 19, 45–56.
#' @export
#' @importFrom Matrix colMeans
#' @importFrom sparseMatrixStats colVars
#' @importFrom utils head tail
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' sfe <- sfe[,sfe$in_tissue]
#' sfe <- logNormCounts(sfe)
#' inds <- order(rowSums(logcounts(sfe)), decreasing = TRUE)[1:50]
#' mat <- logcounts(sfe)[inds,]
#' g <- findVisiumGraph(sfe)
#' out <- multispati_rsp(t(mat), listw = g, nfposi = 10, nfnega = 10)
multispati_rsp <- function(x, listw, nfposi = 30L, nfnega = 30L, scale = TRUE) {
    nfposi <- as.integer(nfposi)
    nfnega <- as.integer(nfnega)
    if (any(c(nfposi, nfnega) < 0L))
        stop("nfposi and nfnega cannot be negative.")
    else if (nfposi == 0L && nfnega == 0L) {
        stop("At least one of nfposi and nfnega must be positive.")
    }
    x <- sweep(x, 2, colMeans(x))
    if (is(x, "dgeMatrix")) x <- as.matrix(x)
    if (scale) {
        # Note that dudi.pca divides by n instead of n-1 when scaling data
        n <- nrow(x)
        x <- sweep(x, 2, sqrt(colVars(x)*(n-1)/n), FUN = "/")
    }
    if (inherits(listw, "Matrix") || is.matrix(listw))
        W <- listw
    else if (is(listw, "listw"))
        W <- listw2sparse(listw)
    else
        stop("listw must be either a listw object or an adjacency matrix.")
    covar <- t(x) %*% (W + t(W)) %*% x / (2*nrow(x))
    if (nfnega == 0L) {
        res <- eigs_sym(covar, k = nfposi, which = "LR")
    } else if (nfposi == 0L) {
        res <- eigs_sym(covar, k = nfnega, which = "SR")
    } else {
        nf <- max(nfposi, nfnega)
        res <- eigs_sym(covar, k = 2*nf, which = "BE")
        inds <- seq_len(2*nf)
        if (nfposi != nfnega) {
            res$values <- c(head(res$values, nfposi), tail(res$values, nfnega))
            res$vectors <- cbind(res$vectors[,head(inds, nfposi)],
                                 res$vectors[,tail(inds, nfnega)])
        }
    }
    loadings <- res$vectors
    out <- x %*% loadings
    colnames(out) <- colnames(loadings) <- paste0("PC", seq_len(ncol(out)))
    rownames(loadings) <- colnames(x)
    attr(out, "rotation") <- loadings
    attr(out, "eig") <- res$values
    out
}

.reorg_multispati <- function(out, x) {
    loadings <- out$vectors
    o <- x %*% loadings
    colnames(o) <- paste0("CS", seq_len(ncol(o)))
    attr(o, "rotation") <- loadings
    attr(o, "eig") <- out$values
    o
}

multispati <- SFEMethod(
    c(name = "multispati", title = "MULTISPATI PCA", package = "Voyager",
      variate = "multi", scope = "global", default_attr = NA),
    fun = multispati_rsp,
    reorganize_fun = function(out) out,
    joint = TRUE
)

# Multivariate version of Local Geary's C=============
#' @include res2df.R

.localC_multi_fun <- function(perm = FALSE) {
    function(x, listw, ..., zero.policy) {
        x <- as.matrix(x)
        fun <- if (perm) spdep::localC_perm else spdep::localC
        fun(x, listw = listw, zero.policy = zero.policy, ...)
    }
}

localC_multi <- SFEMethod(
    c(name = "localC_multi", title = "Multivariate local Geary's C",
      package = "spdep", variate = "multi",
      default_attr = "localC"),
    fun = .localC_multi_fun(FALSE),
    reorganize_fun = .to_df_identity,
    dest = "colData"
)

localC_perm_multi <- SFEMethod(
    c(name = "localC_perm_multi", title = "Multivariate local Geary's C permutation testing",
      package = "spdep", variate = "multi",
      default_attr = "localC"),
    fun = .localC_multi_fun(TRUE),
    reorganize_fun = .localCpermmulti2df,
    dest = "colData"
)
