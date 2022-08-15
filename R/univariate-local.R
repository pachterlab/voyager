# 4. localmoran: Output dataframe can be stored in colData
# 7. Getis-Ord Gi and Gi*
# 11. LOSH
# Where to store the results?
# I suppose colData for gene expression, new column for colData, rowData, and *Geometries

#' Calculate univariate local spatial autocorrelation
#'
#' For global metrics, the entire dataset gets one single value. For local
#' metrics, each location get a value.
#'
#' @inheritParams calculateMoransI
#' @param ... Other arguments passed to the underlying functions: For
#'   \code{calculateLocalI}, \code{\link{localmoran}} or
#'   \code{\link{localmoran_perm}} when \code{perm = TRUE}.
#' @name calculateLocalI
#' @aliases calculateLocalC calculateLocalG calculateLOSH
#' @importFrom spdep localmoran localmoran_perm localC localC_perm localG
#'   localG_perm LOSH LOSH.mc
NULL

.calc_local_fun <- function(fun, fun_perm) {
  function(x, listw, BPPARAM = SerialParam(),
           zero.policy = NULL, perm = FALSE,
           ...) {
    fun <- if (perm) fun_perm else fun
    .calc_univar_autocorr(x, listw, fun = fun, BPPARAM = BPPARAM,
                          returnDF = FALSE, ...)
  }
}

#' @rdname calculateLocalI
#' @export
setMethod("calculateLocalI", "ANY", .calc_local_fun(localmoran, localmoran_perm))

#' @rdname calculateLocalI
#' @export
setMethod("calculateLocalI", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateLocalI))

#' @rdname calcualteLocalI
#' @export
setMethod("calculateLocalC", "ANY", .calc_local_fun(localC, localC_perm))

#' @rdname calculateLocalI
#' @export
setMethod("calculateLocalC", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateLocalC))

#' @rdname calcualteLocalI
#' @export
setMethod("calculateLocalG", "ANY", .calc_local_fun(localG, localG_perm))

#' @rdname calculateLocalI
#' @export
setMethod("calculateLocalG", "SpatialFeatureExperiment",
          .calc_univar_sfe_fun(calculateLocalG))

#' @rdname calculateLocalI
#' @export
calculateLocalGstar <- .calc_univar_sfe_fun(calculateLocalG, returnDF = FALSE,
                                            include_self = TRUE)

#' @rdname calculateLocalI
#' @export
colDataLocalI <- .coldata_univar_fun(calculateLocalI, name = "LocalI")

#' @rdname calculateLocalI
#' @export
colGeometryLocalI <- .colgeom_univar_fun(calculateLocalI, name = "LocalI")

#' @rdname calculateLocalI
#' @export
annotGeometryLocalI <- .annotgeom_univar_fun(calculateLocalI, name = "LocalI")

#' @rdname calculateLocalI
#' @export
runLocalI <- function() {

}
