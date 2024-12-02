#' @include SFEMethod-class.R
#' @include res2df.R

# Global bivariate===================
.scale_n <- function(X) {
    # Divide by n rather than n-1
    X <- sweep(X, 1, rowMeans(X))
    ss <- sqrt(rowSums(X^2))
    sweep(X, 1, ss, "/")
}
#' @importFrom Matrix rowSums rowMeans
#' @importFrom DelayedArray sweep
.lee_mat <- function(x, y = NULL, listw, zero.policy = TRUE, ...) {
    # X has genes in rows
    if (inherits(listw, "listw"))
        W <- listw2sparse(listw)
    else W <- listw
    x <- .scale_n(x)
    if (!is.null(y)) {
        y <- .scale_n(y)
    } else y <- x
    n <- ncol(x) # dimension of y is checked in calculateBivariate
    out <- x %*% (t(W) %*% W) %*% t(y)/sum(rowSums(W)^2) * n
    if (all(dim(out) == 1L)) out <- out[1,1]
    out
}

lee <- SFEMethod(name = "lee", fun = .lee_mat, title = "Lee's bivariate statistic",
                 reorganize_fun = function(out, name, ...) out,
                 package = "Voyager", variate = "bi", scope = "global",
                 use_matrix = TRUE)

lee.mc <- SFEMethod(name = "lee.mc", fun = spdep::lee.mc,
                    title = "Lee's bivariate static with permutation testing",
                    reorganize_fun = .mcsim2df, package = "spdep",
                    variate = "bi", scope = "global", use_matrix = FALSE,
                    args_not_check = "nsim")

lee.test <- SFEMethod(name = "lee.test", fun = spdep::lee.test,
                      title = "Lee's L test", reorganize_fun = .htest2df,
                      package = "spdep", variate = "bi", scope = "global",
                      use_matrix = FALSE)

# Local bivariate================
locallee <- SFEMethod(name = "locallee",
                      fun = function(x, y, listw, zero.policy = TRUE) {
                          spdep::lee(x, y, listw, n = length(x))
                      },
                      title = "Local Lee's bivariate statistic",
                      reorganize_fun = .lee2df, package = "spdep",
                      variate = "bi", scope = "local", use_matrix = FALSE)

localmoran_bv <- SFEMethod(name = "localmoran_bv",
                           fun = function(x, y, listw, ..., zero.policy = TRUE)
                               spdep::localmoran_bv(x, y, listw, ...),
                           title = "Local bivariate Moran's I",
                           reorganize_fun = .LOSHmc2df, package = "spdep",
                           default_attr = "Ibvi",
                           variate = "bi", scope = "local", use_matrix = FALSE)
