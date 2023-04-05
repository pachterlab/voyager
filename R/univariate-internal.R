# Internal function for univariate metrics
.calc_univar <- function(x, fun, BPPARAM, ...) {
    if (is(x, "DFrame") || is.data.frame(x)) {
        if (is(x, "sf")) x <- st_drop_geometry(x)
        x <- t(as.matrix(x))
        if (anyNA(x)) {
            stop("Only numeric columns without NA (within the sample_id) can be used.")
        }
    }
    out <- bplapply(seq_len(nrow(x)), function(i) {
        fun(x[i, ], ...)
    }, BPPARAM = BPPARAM)
    names(out) <- rownames(x)
    return(out)
}

#' @importFrom stats as.formula terms
.get_coords_df <- function(x, df, sample_id, exprs_values,
                           swap_rownames, ...) {
    if (!is(df, "sf") || st_geometry_type(df, by_geometry = FALSE) != "POINT") {
        if (is(df, "sf")) df <- st_drop_geometry(df)
        # Can't use list columns as regressors
        inds_keep <- vapply(df, is.atomic, FUN.VALUE = logical(1))
        df <- df[,inds_keep]
        if (is(df, "DataFrame")) df <- as.data.frame(df)
        colnames_use <- spatialCoordsNames(x)
        geo <- df2sf(spatialCoords(x)[x$sample_id == sample_id,], colnames_use)
        oth_names <- setdiff(names(df), colnames_use)
        if (length(oth_names)) {
            df <- cbind(df[,oth_names], geo)
        } else df <- geo
    }
    dots <- list(...)
    if ("formula" %in% names(dots)) {
        rgs <- labels(terms(dots[["formula"]]))
        if (length(rgs) && any(!rgs %in% names(df))) {
            rgs <- setdiff(rgs, names(df))
            values <- .get_feature_values(x, rgs, sample_id = sample_id,
                                          exprs_values = exprs_values,
                                          show_symbol = !is.null(swap_rownames),
                                          swap_rownames = swap_rownames)
            df <- cbind(df, values)
        }
    }
    st_as_sf(df)
}

#' @importFrom spdep include.self nb2listw
#' @importFrom sf st_as_sf
#' @importFrom SpatialExperiment spatialCoordsNames
.calc_univar_sfe_fun <- function(type = NULL) {
    fun_use <- function(x, type, features = NULL, colGraphName = 1L,
                        colGeometryName = 1L, sample_id = "all",
                        exprs_values = "logcounts", BPPARAM = SerialParam(),
                        zero.policy = NULL, returnDF = TRUE,
                        include_self = FALSE, p.adjust.method = "BH",
                        swap_rownames = NULL, name = NULL, ...) {
        # Am I sure that I want to use logcounts as the default?
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        out <- lapply(sample_id, function(s) {
            features <- .check_features(x, features, swap_rownames = swap_rownames)[["assay"]]
            mat <- assay(x, exprs_values)[features, colData(x)$sample_id == s, drop = FALSE]
            if (use_graph(type)) {
                listw_use <- colGraph(x, type = colGraphName, sample_id = s)
                if (include_self) {
                    nb2 <- include.self(listw_use$neighbours)
                    listw_use <- nb2listw(nb2)
                }
                o <- calculateUnivariate(mat, listw = listw_use,
                                         type = type,
                                         BPPARAM = BPPARAM,
                                         zero.policy = zero.policy,
                                         returnDF = returnDF, p.adjust.method = p.adjust.method,
                                         name = name, ...
                )
            } else {
                cg <- colGeometry(x, colGeometryName, sample_id = s)
                cg <- .get_coords_df(x, cg, s, exprs_values, swap_rownames, ...)
                o <- calculateUnivariate(mat, coords_df = cg, type = type, BPPARAM = BPPARAM,
                                         returnDF = returnDF, p.adjust.method = p.adjust.method,
                                         name = name, ...)
            }
            o
        })
        names(out) <- sample_id
        if (length(sample_id) == 1L) out <- out[[1]]
        out
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features = NULL, colGraphName = 1L, colGeometryName = 1L,
                 sample_id = "all",
                 exprs_values = "logcounts", BPPARAM = SerialParam(),
                 zero.policy = NULL, returnDF = TRUE,
                 include_self = FALSE, p.adjust.method = "BH",
                 swap_rownames = NULL, name = NULL, ...) {
            fun_use(
                x, type, features, colGraphName, colGeometryName, sample_id,
                exprs_values, BPPARAM, zero.policy, returnDF,
                include_self, p.adjust.method, swap_rownames, name, ...
            )
        }
    }
}

# Function to add univariate local results to SFE object

# What should the output look like
# moran: list, no class
# moran.mc: list, class htest and mc.sim
# moran_bv (spdep devel): list
# EBImoran.mc: list, class htest and mc.sim
# geary: list, no class
# lee: list, no class. second element is local lee. Stylistic differences.
# geary.mc: list, htest and mc.sim
# sp.mantel.mc: list, htest and mc.sim
# jointcount.mc: jclist, lists of htest and mc.sim
# lee.mc: list, htest and mc.sim
# globalG: htest
# sp.correlogram: list, class spcor
# moran.plot: data frame
# localmoran: class localmoran, matrix, with attributes on quadrant. Maybe I'll save them as separate types.
# Or convert them into data frames and cbind
# localmoran_perm: class localmoran, matrix with attributes on quadrant
# localmoran_bv (spdep devel): data frame
# localC: vector without attributes. Why is it so different from localmoran?
# localC_perm: class localC, Just like localmoran_perm. Different authors?
# localG: vector, class localG
# localG_perm: class localG, vector with attributes which includes a matrix.
# I think I'll add that vector to that matrix.
# LOSH: matrix, no class
# LOSH.mc: matrix, class LOSH and mc.sim but very different from moran.mc results,
# with simulated p-values
# skater: list, class skater
# Also some stuff from GWmodel
# gwss: list, class lss, SDF sp component for the results. For each type, should be vector.
# It's not modular. I don't like it.

# .get_fun and .graph_fun have arguments x, type, and sample_id
# .get_fun returns a data frame or matrix with column names that correspond to the features
# .set_fd_fun has arguments x, which, sample_id, name, features, res, params, and
# returns an sfe object. The geometry argument is optional.

.make_univar_fun <- function(.get_fun, .graph_fun, .set_fd_fun, type = NULL,
                             colData = FALSE) {
    function(x, type, features, which = NULL, graphName = 1L,
             sample_id = "all",
             BPPARAM = SerialParam(), zero.policy = NULL,
             include_self = FALSE, p.adjust.method = "BH",
             name = NULL, colGeometryName = NULL,
             annotGeometryName = NULL, reducedDimName = NULL, ...) {
        if (is.character(type)) type <- get(type, mode = "S4")
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        if (is.null(name)) name <- info(type, "name")
        other_args <- list(...)
        # But what if different parameters were used to make the graph for
        # different samples? But why would that be a good idea?
        if (use_graph(type))
            g <- .graph_fun(x, type = graphName, sample_id = sample_id[1])
        else g <- NULL
        params <- c(info(type, c("name", "package")),
                    list(version = packageVersion(info(type, "package")),
                         zero.policy = zero.policy, include_self = include_self,
                         p.adjust.method = p.adjust.method,
                         graph_params = attr(g, "method")), other_args)
        local <- is_local(type)
        if (!local) params$p.adjust.method <- NULL
        old_params <- getParams(x, name, colData = colData, local = local,
                                colGeometryName = colGeometryName,
                                annotGeometryName = annotGeometryName,
                                reducedDimName = reducedDimName)
        .check_old_params(params, old_params, name, args_not_check(type))

        for (s in sample_id) {
            df <- .get_fun(x, which, sample_id = s)
            if (use_graph(type)) {
                listw_use <- .graph_fun(x, type = graphName, sample_id = s)
                if (include_self) {
                    nb2 <- include.self(listw_use$neighbours)
                    listw_use <- nb2listw(nb2)
                }
                res <- calculateUnivariate(df[, features, drop = FALSE],
                                           listw = listw_use, type = type, BPPARAM = BPPARAM,
                                           zero.policy = zero.policy, returnDF = TRUE,
                                           p.adjust.method = p.adjust.method, name = name, ...
                )
            } else {
                geo <- .get_coords_df(x, df, s, exprs_values = "logcounts",
                                      swap_rownames = NULL, ...)
                res <- calculateUnivariate(df[, features, drop = FALSE],
                                           coords_df = geo, type = type,
                                           BPPARAM = BPPARAM, returnDF = TRUE,
                                           p.adjust.method = p.adjust.method,
                                           name = name, ...)
            }
            if (local) {
                x <- .add_localResults_info(x, sample_id = s,
                                            name = name, features = features,
                                            res = res, params = params,
                                            colGeometryName = colGeometryName,
                                            annotGeometryName = annotGeometryName)
            } else {
                x <- .set_fd_fun(x, which = which, sample_id = s, name = name,
                                 features = features, res = res,
                                 params = params)
            }
        }
        x
    }
}

.coldata_univar_fun <- function(type = NULL) {
    .get_fun <- function(x, type, sample_id) {
        colData(x)[colData(x)$sample_id == sample_id, , drop = FALSE]
    }
    .set_fd_fun <- function(x, which, sample_id, name, features, res, params)
        .add_fd_dimData(x, MARGIN = 2, res = res, features = features,
                        sample_id = sample_id, name = name, params = params)
    fun <- .make_univar_fun(.get_fun, .graph_fun = colGraph,
                            .set_fd_fun = .set_fd_fun,
                            type = type, colData = TRUE)
    if (is.null(type)) {
        function(x, type, features, colGraphName = 1L, sample_id = "all",
                 BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type, features, which = NULL, graphName = colGraphName,
                sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, ...)
        }
    } else {
        function(x, features, colGraphName = 1L, sample_id = "all",
                 BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type = type, features = features, which = NULL,
                graphName = colGraphName, sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, ...)
        }
    }
}

.colgeom_univar_fun <- function(type = NULL) {
    .set_fd_fun <- function(x, which, sample_id, name, features, res, params) {
        colGeometry(x, which, sample_id = "all") <-
            .add_fd(x, df = colGeometry(x, which, sample_id = "all"),
                    sample_id = sample_id, name = name, features = features,
                    res = res, params = params)
        x
    }
    fun <- .make_univar_fun(.get_fun = colGeometry, .graph_fun = colGraph,
                            .set_fd_fun = .set_fd_fun,
                            type = type)
    if (is.null(type)) {
        function(x, type, features, colGeometryName = 1L, colGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type, features, which = colGeometryName, graphName = colGraphName,
                sample_id = sample_id, BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, colGeometryName = colGeometryName, ...)
        }
    } else {
        function(x, features, colGeometryName = 1L, colGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type = type, features = features, which = colGeometryName,
                graphName = colGraphName, sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, colGeometryName = colGeometryName, ...)
        }
    }
}

.annotgeom_univar_fun <- function(type = NULL) {
    .get_fun <- function(x, type, sample_id) {
        ag <- annotGeometry(x, type, sample_id)
        ag <- .rm_empty_geometries(ag, MARGIN = 3)
        ag
    }
    .set_fd_fun <- function(x, which, sample_id, name, features, res, params) {
        annotGeometry(x, which, sample_id = "all") <-
            .add_fd(
                x, df = annotGeometry(x, which, sample_id = "all"),
                sample_id = sample_id, name = name, features = features,
                res = res, params = params
            )
        x
    }
    fun <- .make_univar_fun(.get_fun = .get_fun, .graph_fun = annotGraph,
                            .set_fd_fun = .set_fd_fun,
                            type = type)
    if (is.null(type)) {
        function(x, type, features, annotGeometryName = 1L, annotGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type, features, which = annotGeometryName, graphName = annotGraphName,
                sample_id = sample_id, BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, annotGeometryName = annotGeometryName, ...)
        }
    } else {
        function(x, features, annotGeometryName = 1L, annotGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            fun(x, type = type, features = features, which = annotGeometryName,
                graphName = annotGraphName, sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, annotGeometryName = annotGeometryName, ...)
        }
    }
}

.reddim_univar_fun <- function(type = NULL) {
    .get_fun <- function(x, type, sample_id) {
        as.data.frame(reducedDim(x, type)[colData(x)$sample_id == sample_id, , drop = FALSE])
    }
    .set_fd_fun <- function(x, which, sample_id, name, features, res, params)
        .add_fd_reddim(x, dimred = which, res = res, features = features,
                       sample_id = sample_id, name = name, params = params)
    fun <- .make_univar_fun(.get_fun, .graph_fun = colGraph,
                            .set_fd_fun = .set_fd_fun,
                            type = type)
    if (is.null(type)) {
        function(x, type, dimred = 1L, components = 1L, colGraphName = 1L,
                 sample_id = "all",
                 BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            x <- .add_reddim_colnames(x, dimred)
            features <- colnames(reducedDim(x, dimred))[components]
            fun(x, type, features = features, which = dimred, graphName = colGraphName,
                sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, reducedDimName = dimred, ...)
        }
    } else {
        function(x, dimred = 1L, components = 1L, colGraphName = 1L, sample_id = "all",
                 BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH",
                 name = NULL, ...) {
            x <- .add_reddim_colnames(x, dimred)
            features <- colnames(reducedDim(x, dimred))[components]
            fun(x, type = type, features = features, which = dimred,
                graphName = colGraphName, sample_id = sample_id,
                BPPARAM = BPPARAM, zero.policy = zero.policy,
                include_self = include_self, p.adjust.method = p.adjust.method,
                name = name, reducedDimName = dimred, ...)
        }
    }
}

.sfe_univar_fun <- function(type = NULL) {
    fun_use <- function(x, type, features = NULL, colGraphName = 1L,
                        colGeometryName = 1L, sample_id = "all",
                        exprs_values = "logcounts", BPPARAM = SerialParam(),
                        swap_rownames = NULL,
                        zero.policy = NULL, include_self = FALSE,
                        p.adjust.method = "BH", name = NULL, ...) {
        if (is.character(type)) type <- get(type, mode = "S4")
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        if (is.null(name)) name <- info(type, "name")
        other_args <- list(...)
        if (use_graph(type))
            g <- colGraph(x, type = colGraphName, sample_id = sample_id[1])
        else g <- NULL
        params <- c(info(type, c("name", "package")),
                    list(version = packageVersion(info(type, "package")),
                         zero.policy = zero.policy, include_self = include_self,
                         p.adjust.method = p.adjust.method,
                         graph_params = attr(g, "method")), other_args)
        local <- is_local(type)
        if (!local) params$p.adjust.method <- NULL
        old_params <- getParams(x, name, local = local)
        .check_old_params(params, old_params, name, args_not_check(type))

        if (is.null(features)) features <- rownames(x)
        # Requires devel version of SFE
        features <- .symbol2id(x, features, swap_rownames)
        for (s in sample_id) {
            out <- calculateUnivariate(x, type, features, colGraphName, colGeometryName, s,
                                       exprs_values, BPPARAM, zero.policy,
                                       returnDF = TRUE,
                                       include_self = include_self, p.adjust.method = p.adjust.method,
                                       ...
            )
            if (local) {
                x <- .add_localResults_info(x, sample_id = s,
                                            name = name, features = features,
                                            res = out, params = params)
            } else {
                out <- .add_name_sample_id(out, s)
                rowData(x)[features, names(out)] <- out
                metadata(rowData(x))$params[[name]] <- params
            }
        }
        x
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features = NULL, colGraphName = 1L, colGeometryName = 1L,
                 sample_id = "all",
                 exprs_values = "logcounts", BPPARAM = SerialParam(),
                 swap_rownames = NULL,
                 zero.policy = NULL, include_self = FALSE,
                 p.adjust.method = "BH", name = NULL, ...) {
            fun_use(
                x, type, features, colGraphName, colGeometryName, sample_id,
                exprs_values, BPPARAM, swap_rownames, zero.policy, include_self,
                p.adjust.method, name, ...
            )
        }
    }
}
