#' @export
setGeneric("calculateUnivariate", function(x, type, ...) standardGeneric("calculateUnivariate"))

#' @export
setGeneric("calculateMoransI", function(x, ...) standardGeneric("calculateMoransI"))

# @export
# setGeneric("calculateBivariate", function(x, ...) standardGeneric("calculateBivariate"))

# @export
# setGeneric("calculateMultivariate", function(x, ...) standardGeneric("calculateMultivariate"))

#' @export
setGeneric("info", function(x, type) standardGeneric("info"))

#' @export
setGeneric("fun", function(x) standardGeneric("fun"))

#' @export
setGeneric("to_df_fun", function(x) standardGeneric("to_df_fun"))

#' @export
setGeneric("is_local", function(x) standardGeneric("is_local"))
