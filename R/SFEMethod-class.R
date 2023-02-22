#' SFEMethod class
#'
#' This S4 class is used to wrap spatial analysis methods, taking inspiration
#' from the \code{caret} and \code{tidymodels} packages.
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
#' @slot fun The function implementing the method. For univariate methods, the
#'   first two arguments must be \code{x} for a vector, and \code{listw} for a
#'   \code{listw} object specifying the spatial neighborhood graph, and
#'   optionally \code{...}. If the original function implementing the method in
#'   the package has different argument names or orders, write a thin wrapper to
#'   rearrange and/or rename the arguments.
#' @slot to_df_fun Function to convert output from \code{fun} into a format to
#'   store in the SFE object. For univariate global methods, different fields of
#'   the result should be columns of a data frame with one row so results for
#'   multiple features will be a data frame. The arguments should be \code{out}
#'   for a list of raw output, each element of which is output for one feature
#'   and \code{name} to rename the primary field if a more informative name is
#'   needed, and \code{...} for other arguments specific to methods. For
#'   univariate local methods, the output should be a data frame or matrix whose
#'   rows match the columns of the gene count matrix. The arguments should be
#'   \code{out}, \code{nb} for a neighborhood list used for multiple testing
#'   correction, and \code{p.adjust.method} for a method to correct for multiple
#'   testing as in \code{\link{p.adjust}}, and \code{...}.
#' @name SFEMethod
#' @aliases SFEMethod-class
setClass("SFEMethod", slots = c(
    info = "character",
    fun = "function",
    to_df_fun = "function"
))

.valid_SFEMethod <- function(object) {
    outs <- list()
    nms <- c("name", "variate", "scope", "package", "title", "default_attr")
    if (!setequal(names(object@info), nms))
        outs <- c(outs, paste("Slot `info` must have names",
                              paste(nms, collapse = ", ")))
    variates <- c("uni", "bi", "multi")
    if (!object@info["variate"] %in% variates) {
        outs <- c(outs, paste("Field 'variate' in slot `info` must be one of",
                              paste(variates, collapse = ", ")))
    }
    scopes <- c("global", "local")
    if (!object@info["scope"] %in% scopes) {
        outs <- c(outs, paste("Field 'scope' in slot `info` must be one of",
                              paste(scopes, collapse = " and ")))
    }
    if (object@info["variate"] == "univariate") {
        if (names(formals(object@fun)) != c("x", "listw", "..."))
            outs <- c(outs, "The first two arguments of slot `fun` must be 'x' and 'listw'")
        if (object@info["scope"] == "global") {
            if (names(formals(object@to_df_fun)) != c("out", "name", "..."))
                outs <- c(outs, "Slot `to_df_fun` must have arguments 'out', 'name', and '...'")
        } else {
            if (names(formals(object@to_df_fun)) != c("out", "nb", "p.adjust.method", "..."))
                outs <- c(outs, "Slot `to_df_fun` must have arguments 'out', 'nb', 'p.adjust.method', and '...'")
        }
    }
    if (!length(outs)) return(TRUE)
    outs
}

setValidity("SFEMethod", .valid_SFEMethod)

# I don't think I need setters yet, so only getters for now.

#' @param info See slot documentation
#' @param fun See slot documentation
#' @param to_df_fun See slot documentation
#' @param x A \code{SFEMethod} object
#' @param type One of the names of the \code{info} slot, see slot documentation.
#' @return The constructor returns a \code{SFEMethod} object. The getters return
#' the content of the corresponding slots.
#' @export
#' @rdname SFEMethod
SFEMethod <- function(info, fun, to_df_fun) {
    new("SFEMethod", info = info, fun = fun, to_df_fun = to_df_fun)
}

#' @export
#' @rdname SFEMethod
setMethod("info", "SFEMethod", function(x, type) x@info[type])

#' @export
#' @rdname SFEMethod
setMethod("is_local", "SFEMethod", function(x) info(x, "scope") == "local")

#' @export
#' @rdname SFEMethod
setMethod("fun", "SFEMethod", function(x) x@fun)

#' @export
#' @rdname SFEMethod
setMethod("to_df_fun", "SFEMethod", function(x) x@to_df_fun)
