.check_old_params <- function(params, old_params, name, not_check) {
    if (!identical(not_check, NA_character_))
        params <- params[!names(params) %in% not_check]
    # Maybe temporary fix, but I expect graphs of the same name for different
    # samples to come from the same parameters.
    params$graph_params$args$sample_id <- NULL
    old_params$graph_params$args$sample_id <- NULL
    if (length(old_params)) {
        old_params <- old_params[names(params)]
        is_matched <- vapply(seq_along(params), function(i) {
            isTRUE(all.equal(params[[i]], old_params[[i]]))
        }, FUN.VALUE = logical(1L))
        nms <- c("version", "package", "graph_params")
        inds <- names(params) %in% nms
        if (any(!is_matched[!inds]) && all(is_matched[inds])) {
            stop("New results were computed with different parameters ",
                 "from existing results in localResult ", name,
                 "; please use a different name for different parameters.")
        }
        if (!is_matched[names(params) == "package"]) {
            message("New results were computed with package ", params[["package"]],
                    " while existing results used ", old_params[["package"]],
                    "; please verify consistency between packages.")
        } else if (!is_matched[names(params) == "version"]) {
            message("New results were computed with version ",
                    params[["version"]], " of package ",
                    params[["package"]], ", while existing results used ",
                    "version ", old_params[["version"]], " in localResult ",
                    name, "; please verify consistency between versions.")
        }
        if (!is_matched[names(params) == "graph_params"]) {
            # Not throwing errors here since implementations in different packages
            # can give consistent results. Just to notify users to double check.
            if (length(old_params$graph_params)) {
                is_matched_graph <- vapply(
                    seq_along(params$graph_params),
                    function(i) {
                        isTRUE(all.equal(params$graph_params[[i]],
                                         old_params$graph_params[[i]]))
                    }, FUN.VALUE = logical(1L))
            }
            if (!is_matched_graph[[2]]) { # package
                if (params$graph_params$package[[1]] !=
                    old_params$graph_params$package[[1]]) {
                    message("New results used package ", params$graph_params$package[[1]],
                            " to compute spatial neighborhood graph while existing results used ",
                            old_params$graph_params$package[[1]])
                } else if (length(old_params$graph_params$package) > 1L) {
                    message("New results used version ",
                            as.character(params$graph_params$package[[2]]),
                            " of package ", params$graph_params$package[[1]],
                            " to compute spatial neighborhood graph while existing results used ",
                            as.character(old_params$graph_params$package[[2]]))
                }
            } else if (!is_matched_graph[[1]]) { # FUN
                message("New results used a spatial neighborhood graph computed with function ",
                        params$graph_params$FUN, " while existing results used ",
                        old_params$graph_params$FUN)
            } else {
                message("New results come from a spatial neighborhood graph that ",
                        "was computed with different parameters from that used in existing results.")
            }
        }
    }
}

#' @importFrom S4Vectors make_zero_col_DFrame
.initialize_featureData <- function(df) {
    if (is.null(attr(df, "featureData"))) {
        fd <- make_zero_col_DFrame(nrow = ncol(df))
        rownames(fd) <- colnames(df)
        attr(df, "featureData") <- fd
    }
    df
}

.add_fd <- function(x, df, sample_id, name, features, res, params) {
    res <- .add_name_sample_id(res, sample_id)
    df <- .initialize_featureData(df)
    fd <- attr(df, "featureData")
    fd[features, names(res)] <- res
    attr(df, "featureData") <- fd
    attr(df, "params")[[name]] <- params
    df
}

#' @importFrom S4Vectors metadata metadata<-
.initialize_fd_dimData <- function(x, MARGIN) {
    fd_name <- "featureData"
    dimData <- switch(MARGIN, rowData, colData)
    `dimData<-` <- switch(MARGIN, `rowData<-`, `colData<-`)
    if (is.null(metadata(dimData(x))[[fd_name]])) {
        rownames_use <- names(dimData(x))
        fd <- make_zero_col_DFrame(nrow = length(rownames_use))
        rownames(fd) <- rownames_use
        metadata(dimData(x))[[fd_name]] <- fd
    }
    x
}

.add_fd_dimData <- function(x, MARGIN, sample_id, name, features, res, params) {
    res <- .add_name_sample_id(res, sample_id)
    x <- .initialize_fd_dimData(x, MARGIN)
    fd_name <- "featureData"
    dimData <- switch(MARGIN, rowData, colData)
    `dimData<-` <- switch(MARGIN, `rowData<-`, `colData<-`)
    fd <- metadata(dimData(x))[[fd_name]]
    fd[features, names(res)] <- res
    metadata(dimData(x))[[fd_name]] <- fd
    metadata(dimData(x))$params[[name]] <- params
    x
}

.add_localResults_info <- function(x, sample_id, name, features, res, params,
                                   colGeometryName = NULL,
                                   annotGeometryName = NULL) {
    localResults(x, sample_id, name, features,
                 colGeometryName = colGeometryName,
                 annotGeometryName = annotGeometryName) <- res
    if (is.null(colGeometryName)) {
        if (is.null(annotGeometryName))
            metadata(int_colData(x)$localResults)$params[[name]] <- params
        else {
            attr(annotGeometry(x, annotGeometryName, "all")$localResults, "params")[[name]] <- params
        }
    } else {
        attr(colGeometry(x, colGeometryName, "all")$localResults, "params")[[name]] <- params
    }
    x
}

#' Get metadata of colData, rowData, and geometries
#'
#' Results of spatial analyses on columns in \code{colData}, \code{rowData}, and
#' geometries are stored in their metadata, which can be accessed by the
#' \code{\link{metadata}} function. The \code{colFeaturedata} function allows
#' the users to more directly access these results.
#'
#' @param sfe An SFE object.
#' @param type Which geometry, can be name (character) or index (integer)
#' @param MARGIN Integer, 1 means rowGeometry, 2 means colGeometry, and 3 means
#'   annotGeometry. Defaults to 2, colGeometry.
#' @return A \code{DataFrame}.
#' @seealso getParams
#' @export
#' @name colFeatureData
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Moran's I for colData
#' sfe <- colDataMoransI(sfe, "nCounts")
#' colFeatureData(sfe)
colFeatureData <- function(sfe) {
    metadata(colData(sfe))$featureData
}

#' @rdname colFeatureData
#' @export
rowFeatureData <- function(sfe) {
    metadata(rowData(sfe))$featureData
}

#' @rdname colFeatureData
#' @export
geometryFeatureData <- function(sfe, type, MARGIN = 2L) {
    geo_fun <- switch (MARGIN, rowGeometry, colGeometry, annotGeometry)
    df <- geo_fun(sfe, type, sample_id = "all")
    attr(df, "featureData")
}

#' Get parameters used in spatial methods
#'
#' The \code{getParams} function allows users to access the parameters used to
#' compute the results that may be stored in \code{\link{colFeatureData}}.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param name Name used to store the results.
#' @param local Logical, whether the results of interest come from a local
#'   spatial method.
#' @param colData Logical, whether the results were computed for a column of
#'   \code{colData(sfe)}.
#' @param colGeometryName To get results for a \code{colGeometry}.
#' @param annotGeometryName To get results for an \code{annotGeometry};
#'   \code{colGeometry} has precedence so this argument is ignored if
#'   \code{colGeometryName} is specified.
#' @return A named list showing the parameters
#' @rdname colFeatureData
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' sfe <- colDataMoransI(sfe, "nCounts")
#' getParams(sfe, "moran", colData = TRUE)
getParams <- function(sfe, name, local = FALSE, colData = FALSE,
                      colGeometryName = NULL, annotGeometryName = NULL) {
    if (local) {
        if (is.null(colGeometryName)) {
            if (is.null(annotGeometryName)) {
                lr <- int_colData(sfe)$localResults
                if (is.null(lr)) return(NULL)
                return(metadata(lr)$params[[name]])
            } else {
                ag <- annotGeometry(sfe, annotGeometryName, "all")
                lr <- ag$localResults
                return(attr(lr, "params")[[name]])
            }
        } else {
            cg <- colGeometry(sfe, colGeometryName, "all")
            lr <- cg$localResults
            return(attr(lr, "params")[[name]])
        }
    } else {
        if (colData) {
            metadata(colData(sfe))$params[[name]]
        } else if (is.null(colGeometryName)) {
            if (is.null(annotGeometryName)) {
                metadata(rowData(sfe))$params[[name]]
            } else {
                attr(annotGeometry(sfe, annotGeometryName, sample_id = "all"), "params")[[name]]
            }
        } else {
            attr(colGeometry(sfe, colGeometryName, sample_id = "all"), "params")[[name]]
        }
    }
}
