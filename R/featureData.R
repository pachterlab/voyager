#' @importFrom S4Vectors mcols metadata metadata<- mcols<-
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
        if ("graph_params" %in% names(params) &&
            !is_matched[names(params) == "graph_params"]) {
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

.initDF <- function(m) {
    rownames_use <- colnames(m)
    fd <- make_zero_col_DFrame(nrow = ncol(m))
    rownames(fd) <- rownames_use
    fd
}
.initialize_fd_reddim <- function(x, dimred) {
    if (is.null(attr(reducedDim(x, dimred), "featureData"))) {
        attr(reducedDim(x, dimred), "featureData") <- .initDF(reducedDim(x, dimred))
    }
    x
}

.add_fd_dimData <- function(x, MARGIN, sample_id, name, features, res, params) {
    res <- .add_name_sample_id(res, sample_id)
    dimData <- switch(MARGIN, rowData, colData)
    `dimData<-` <- switch(MARGIN, `rowData<-`, `colData<-`)
    if (is.null(mcols(dimData(x)))) fd <- .initDF(dimData(x))
    else fd <- mcols(dimData(x)) 
    fd[features, names(res)] <- res
    mcols(dimData(x)) <- fd
    metadata(dimData(x))$params[[name]] <- params
    x
}

.add_fd_reddim <- function(x, dimred, sample_id, name, features, res, params) {
    res <- .add_name_sample_id(res, sample_id)
    x <- .initialize_fd_reddim(x, dimred)
    fd <- attr(reducedDim(x, dimred), "featureData")
    fd[features, names(res)] <- res
    attr(reducedDim(x, dimred), "featureData") <- fd
    attr(reducedDim(x, dimred), "params")[[name]] <- params
    x
}

.add_localResults_info <- function(x, sample_id, name, features, res, params,
                                   colGeometryName = NULL,
                                   annotGeometryName = NULL,
                                   reducedDimName = NULL) {
    localResults(x, sample_id, name, features,
                 colGeometryName = colGeometryName,
                 annotGeometryName = annotGeometryName) <- res
    if (is.null(colGeometryName)) {
        if (is.null(annotGeometryName))
            metadata(int_colData(x)$localResults)$params[[name]] <- params
        else if (is.null(reducedDimName)) {
            attr(annotGeometry(x, annotGeometryName, "all")$localResults, "params")[[name]] <- params
        }
    } else {
        attr(colGeometry(x, colGeometryName, "all")$localResults, "params")[[name]] <- params
    }
    x
}
