#' SFEMethod class
#'
#' This S4 class is used to wrap spatial analysis methods, taking inspiration
#' from the \code{caret} and \code{tidymodels} packages.
#'
#' The \code{fun} slot should be specified as such:
#'
#' For all methods, there must be arguments \code{x} for a vector, \code{listw}
#' for a \code{listw} object specifying the spatial neighborhood graph,
#' \code{zero.policy} specifying what to do with cells without neighbors
#' (default NULL, use global option value; if TRUE assign zero to the lagged
#' value of zones without neighbours, if FALSE assign NA), and optionally other
#' method specific arguments and \code{...} to pass to the underlying imported
#' function. If the original function implementing the method in the package has
#' different argument names or orders, write a thin wrapper to rearrange and/or
#' rename the arguments.
#'
#' For univariate methods, the first two arguments must be \code{x} and
#' \code{listw}.
#'
#' For bivariate methods, the first three arguments must be \code{x}, \code{y},
#' and \code{listw}.
#'
#' For multivariate methods, the argument \code{x} is mandatory, for the matrix
#' input. These arguments must be present but can be optional by having
#' defaults: \code{listw} and \code{ncomponents} to set the number of dimentions
#' in the output.
#'
#' The \code{reorganize_fun} slot should be specified as such:
#'
#' For univariate global methods, different fields of the result should be
#' columns of a data frame with one row so results for multiple features will be
#' a data frame. The arguments should be \code{out} for a list of raw output,
#' each element of which is output for one feature and \code{name} to rename the
#' primary field if a more informative name is needed, and \code{...} for other
#' arguments specific to methods.
#'
#' For univariate local methods, the output should be a data frame or matrix
#' whose rows match the columns of the gene count matrix. The arguments should
#' be \code{out}, \code{nb} for a neighborhood list used for multiple testing
#' correction, and \code{p.adjust.method} for a method to correct for multiple
#' testing as in \code{\link{p.adjust}}, and \code{...}.
#'
#' @slot info A named character vector specifying information about the method:
#' \describe{
#' \item{name}{Name of the method, used by user-facing functions to specify the
#' method to use, such as "moran" for Moran's I.}
#' \item{variate}{How many variables this method works with, must be one of
#' "uni" for univariate, "bi" for bivariate, or "multi" for multivariate.}
#' \item{scope}{Either "global", returning one result for tne entire dataset,
#' or "local", returning one result for each spatial location.}
#' \item{package}{Name of the package whose implementation of the method is used
#' here, used to check if the package is installed.}
#' \item{title}{Descriptive title to show when plotting the results.}
#' \item{default_attr}{For local methods that return multiple fields, such as
#' local Moran values and their p-values, the default field to use when
#' plotting.}
#' }
#' @slot fun The function implementing the method. See Details.
#' @slot reorganize_fun Function to convert output from \code{fun} into a format to
#'   store in the SFE object. See Details.
#' @name SFEMethod
#' @aliases SFEMethod-class
setClass("SFEMethod", slots = c(
    info = "character",
    fun = "function",
    reorganize_fun = "function"
))

.valid_SFEMethod <- function(object) {
    outs <- list()
    nms <- c("name", "variate", "scope", "package", "title", "default_attr")
    if (!setequal(names(object@info), nms))
        return(paste("Slot `info` must have names",paste(nms, collapse = ", ")))
    variates <- c("uni", "bi", "multi")
    if (!object@info["variate"] %in% variates) {
        return(paste("Field 'variate' in slot `info` must be one of",
                     paste(variates, collapse = ", ")))
    }
    scopes <- c("global", "local")
    if (!object@info["scope"] %in% scopes) {
        return(paste("Field 'scope' in slot `info` must be one of",
                     paste(scopes, collapse = " and ")))
    }
    if (object@info["variate"] == "uni") {
        fm <- names(formals(object@fun))
        if (!identical(fm[seq_len(2)], c("x", "listw")))
            outs <- c(outs, "The first two arguments of slot `fun` must be 'x' and 'listw'")
        if (!"zero.policy" %in% fm)
            outs <- c(outs, "zero.policy must be an argument of slot `fun`")
        if (object@info["scope"] == "global") {
            if (!identical(names(formals(object@reorganize_fun)), c("out", "name", "...")))
                outs <- c(outs, "Slot `reorganize_fun` must have arguments 'out', 'name', and '...'")
        } else {
            if (!identical(names(formals(object@reorganize_fun)), c("out", "nb", "p.adjust.method", "call")))
                outs <- c(outs, "Slot `reorganize_fun` must have arguments 'out', 'nb', 'p.adjust.method', and 'call'")
        }
    }
    if (!length(outs)) return(TRUE)
    unlist(outs)
}

setValidity("SFEMethod", .valid_SFEMethod)

# I don't think I need setters yet, so only getters for now.

#' @param info See slot documentation
#' @param fun See Details.
#' @param reorganize_fun See Details.
#' @param x A \code{SFEMethod} object
#' @param type One of the names of the \code{info} slot, see slot documentation.
#' @return The constructor returns a \code{SFEMethod} object. The getters return
#' the content of the corresponding slots.
#' @export
#' @rdname SFEMethod
SFEMethod <- function(info, fun, reorganize_fun) {
    new("SFEMethod", info = info, fun = fun, reorganize_fun = reorganize_fun)
}

#' @export
#' @rdname SFEMethod
setMethod("info", "SFEMethod", function(x, type) {
    if (missing(type)) x@info else x@info[type]
})

#' @export
#' @rdname SFEMethod
setMethod("is_local", "SFEMethod", function(x) info(x, "scope") == "local")

#' @export
#' @rdname SFEMethod
setMethod("fun", "SFEMethod", function(x) x@fun)

#' @export
#' @rdname SFEMethod
setMethod("reorganize_fun", "SFEMethod", function(x) x@reorganize_fun)
