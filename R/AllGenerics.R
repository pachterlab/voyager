#' @export
setGeneric("calculateUnivariate", function(x, type, ...) standardGeneric("calculateUnivariate"))

#' @export
setGeneric("calculateMoransI", function(x, ...) standardGeneric("calculateMoransI"))

#' @export
setGeneric("calculateBivariate", function(x, ...) standardGeneric("calculateBivariate"))

#' @export
setGeneric("calculateMultivariate", function(x, type, ...) standardGeneric("calculateMultivariate"))

#' @export
setGeneric("info", function(x, type) standardGeneric("info"))

#' @export
setGeneric("fun", function(x) standardGeneric("fun"))

#' @export
setGeneric("reorganize_fun", function(x) standardGeneric("reorganize_fun"))

#' @export
setGeneric("is_local", function(x) standardGeneric("is_local"))

#' @export
setGeneric("args_not_check", function(x) standardGeneric("args_not_check"))

#' @export
setGeneric("is_joint", function(x) standardGeneric("is_joint"))

#' @export
setGeneric("use_graph", function(x) standardGeneric("use_graph"))

#' @export
setGeneric("use_matrix", function(x) standardGeneric("use_matrix"))
