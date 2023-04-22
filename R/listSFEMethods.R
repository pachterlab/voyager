#' List all spatial methods in Voyager package
#'
#' This package ships with many spatial statistics methods as
#' \code{\link{SFEMethod}} objects. The user can adapt the uniform user
#' interface of this package to other spatial methods by creating new
#' \code{SFEMethod} objects. This function lists the names of all methods within
#' \code{Voyager}, to use for the \code{type} argument in
#' \code{\link{calculateUnivariate}}, \code{\link{calculateBivariate}}, and
#' \code{\link{calculateMultivariate}}.
#'
#' @param variate Uni-, bi-, or multi-variate.
#' @param scope whether it's local or global.
#' @return A data frame with a column for the name and another for a brief
#'   description.
#' @export
#' @examples
#' listSFEMethods("uni", "local")
listSFEMethods <- function(variate = c("uni", "bi", "multi"),
                        scope = c("global", "local")) {
    variate <- match.arg(variate)
    scope <- match.arg(scope)
    if (variate == "multi")
        out <- sfe_methods[sfe_methods$variate == variate,]
    else
        out <- sfe_methods[sfe_methods$variate == variate & sfe_methods$scope == scope,]
    out$variate <- NULL
    out$scope <- NULL
    rownames(out) <- NULL
    out
}
