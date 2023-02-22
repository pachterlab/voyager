# Internal function for univariate metrics
.calc_univar <- function(x, listw, fun, BPPARAM, ...) {
    if (is.null(listw)) {
        stop("The graph specified is absent from the SFE object.")
    }
    if (is.vector(x)) {
        x <- matrix(x, nrow = 1)
    }
    if (is(x, "DFrame") || is.data.frame(x)) {
        if (is(x, "sf")) x <- st_drop_geometry(x)
        x <- t(as.matrix(x))
        if (anyNA(x)) {
            stop("Only numeric columns without NA (within the sample_id) can be used.")
        }
    }
    out <- bplapply(seq_len(nrow(x)), function(i) {
        fun(x[i, ], listw, ...)
    }, BPPARAM = BPPARAM)
    names(out) <- rownames(x)
    return(out)
}

.obscure_arg_defaults <- function(listw, type) {
    nb <- listw$neighbours
    switch(type,
        moran = list(n = length(nb), S0 = Szero(listw)),
        geary = list(
            n = length(nb), n1 = length(nb) - 1,
            S0 = Szero(listw)
        ),
        lee = list(n = length(nb)),
        sp.correlogram = list(method = "I"),
        moran.plot = list(plot = FALSE)
    )
}

#' @importFrom spdep include.self nb2listw
.calc_univar_sfe_fun <- function(type = NULL) {
    fun_use <- function(x, type, features = NULL, colGraphName = 1L,
                        sample_id = "all",
                        exprs_values = "logcounts", BPPARAM = SerialParam(),
                        zero.policy = NULL, returnDF = TRUE,
                        include_self = FALSE, p.adjust.method = "BH",
                        swap_rownames = NULL, ...) {
        # Am I sure that I want to use logcounts as the default?
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        out <- lapply(sample_id, function(s) {
            features <- .check_features(x, features, swap_rownames = swap_rownames)[["assay"]]
            listw_use <- colGraph(x, type = colGraphName, sample_id = s)
            if (include_self) {
                nb2 <- include.self(listw_use$neighbours)
                listw_use <- nb2listw(nb2)
            }
            mat <- assay(x, exprs_values)[features, colData(x)$sample_id == s]
            o <- calculateUnivariate(mat, listw = listw_use,
                type = type,
                BPPARAM = BPPARAM,
                zero.policy = zero.policy,
                returnDF = returnDF, p.adjust.method = p.adjust.method, ...
            )
            o
        })
        names(out) <- sample_id
        if (length(sample_id) == 1L) out <- out[[1]]
        out
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features = NULL, colGraphName = 1L, sample_id = "all",
                 exprs_values = "logcounts", BPPARAM = SerialParam(),
                 zero.policy = NULL, returnDF = TRUE,
                 include_self = FALSE, p.adjust.method = "BH",
                 swap_rownames = NULL, ...) {
            fun_use(
                x, type, features, colGraphName, sample_id,
                exprs_values, BPPARAM, zero.policy, returnDF,
                include_self, p.adjust.method, swap_rownames, ...
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
# There's a lot in ESDA. I no longer feel like writing a separate function for
# each method. Shall I refactor the code so everything univariate global will be
# calculateUnivar or runUnivar and everything univariate local will be
# calculateUnivarLocal or runUnivarLocal, and everything multivariate will be
# runSpatialDimRed? And common methods like Moran's I and Gi* will be special
# cases. So I can rename those now internal functions and export them.
# Maybe bivariate should be its own category separate from multivariate.
# For bivariate, shall I force users to do one pair at a time, or do all pairwise
# combinations of a vector of features, or use a list of pairs (can be matrix)?
# I suppose for the first version, I'll focus on univariate.

# Also need plotting functions for localResults
# And functions to compute the spatial metrics for reducedDims
# And plot reducedDims values in space

.is_local <- function(type) {
    if (type %in% c(
        "localmoran", "localmoran_perm", "localC", "localC_perm",
        "localG", "localG_perm", "LOSH", "LOSH.mc", "LOSH.cs", "gwss",
        "lee", "localmoran_bv", "moran.plot"
    )) {
        TRUE
    } else {
        FALSE
    }
}

.coldata_univar_fun <- function(type = NULL) {
    fun_use <- function(x, type, features, colGraphName = 1L, sample_id = "all",
                        BPPARAM = SerialParam(), zero.policy = NULL,
                        include_self = FALSE, p.adjust.method = "BH", ...) {
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        for (s in sample_id) {
            listw_use <- colGraph(x, type = colGraphName, sample_id = s)
            if (include_self) {
                nb2 <- include.self(listw_use$neighbours)
                listw_use <- nb2listw(nb2)
            }
            res <- calculateUnivariate(colData(x)[colData(x)$sample_id == s,
                                                  features, drop = FALSE],
                listw = listw_use, type = type, BPPARAM = BPPARAM,
                zero.policy = zero.policy, returnDF = TRUE,
                p.adjust.method = p.adjust.method, ...
            )
            local <- .is_local(type)
            if (local) {
                localResults(x, s, type, features) <- res
            } else {
                x <- .add_fd_dimData(x, MARGIN = 2, res, features, s, type, ...)
            }
        }
        x
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features, colGraphName = 1L, sample_id = "all",
                 BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH", ...) {
            fun_use(
                x, type, features, colGraphName, sample_id,
                BPPARAM, zero.policy, include_self, p.adjust.method, ...
            )
        }
    }
}

.colgeom_univar_fun <- function(type = NULL) {
    fun_use <- function(x, type, features, colGeometryName = 1L,
                        colGraphName = 1L, sample_id = "all",
                        BPPARAM = SerialParam(), zero.policy = NULL,
                        include_self = FALSE, p.adjust.method = "BH", ...) {
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        for (s in sample_id) {
            listw_use <- colGraph(x, type = colGraphName, sample_id = s)
            if (include_self) {
                nb2 <- include.self(listw_use$neighbours)
                listw_use <- nb2listw(nb2)
            }
            cg <- colGeometry(x, type = colGeometryName, sample_id = s)
            res <- calculateUnivariate(cg[, features, drop = FALSE],
                                       listw = listw_use,
                type = type, BPPARAM = BPPARAM, zero.policy = zero.policy,
                returnDF = TRUE, p.adjust.method = p.adjust.method, ...
            )
            local <- .is_local(type)
            if (local) {
                localResults(x, s, type, features,
                    colGeometryName = colGeometryName
                ) <- res
            } else {
                colGeometry(x, colGeometryName, sample_id = "all") <-
                    .add_fd(
                        x, colGeometry(x, colGeometryName, sample_id = "all"),
                        res, features, s, type
                    )
            }
        }
        x
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features, colGeometryName = 1L, colGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH", ...) {
            fun_use(
                x, type, features, colGeometryName, colGraphName, sample_id,
                BPPARAM, zero.policy, include_self, p.adjust.method, ...
            )
        }
    }
}

.annotgeom_univar_fun <- function(type = NULL) {
    fun_use <- function(x, type, features, annotGeometryName = 1L,
                        annotGraphName = 1L, sample_id = "all",
                        BPPARAM = SerialParam(), zero.policy = NULL,
                        include_self = FALSE, p.adjust.method = "BH", ...) {
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        for (s in sample_id) {
            listw_use <- annotGraph(x, type = annotGraphName, sample_id = s)
            if (include_self) {
                nb2 <- include.self(listw_use$neighbours)
                listw_use <- nb2listw(nb2)
            }
            ag <- annotGeometry(x, type = annotGeometryName, sample_id = s)
            ag <- .rm_empty_geometries(ag, MARGIN = 3)
            res <- calculateUnivariate(
                ag[, features, drop = FALSE], listw = listw_use,
                type = type, BPPARAM = BPPARAM, zero.policy = zero.policy,
                returnDF = TRUE, p.adjust.method = p.adjust.method, ...
            )
            local <- .is_local(type)
            if (local) {
                localResults(x, s, type, features,
                    annotGeometryName = annotGeometryName
                ) <- res
            } else {
                annotGeometry(x, annotGeometryName, sample_id = "all") <-
                    .add_fd(
                        x, annotGeometry(x, annotGeometryName,
                                         sample_id = "all"),
                        res, features, s, type
                    )
            }
        }
        x
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features, annotGeometryName = 1L, annotGraphName = 1L,
                 sample_id = "all", BPPARAM = SerialParam(), zero.policy = NULL,
                 include_self = FALSE, p.adjust.method = "BH", ...) {
            fun_use(
                x, type, features, annotGeometryName, annotGraphName, sample_id,
                BPPARAM, zero.policy, include_self, p.adjust.method, ...
            )
        }
    }
}

.sfe_univar_fun <- function(type = NULL) {
    fun_use <- function(x, type, features = NULL, colGraphName = 1L, sample_id = "all",
                        exprs_values = "logcounts", BPPARAM = SerialParam(),
                        swap_rownames = NULL,
                        zero.policy = NULL, include_self = FALSE,
                        p.adjust.method = "BH", ...) {
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        if (is.null(features)) features <- rownames(x)
        # Requires devel version of SFE
        features <- .symbol2id(x, features, swap_rownames)
        for (s in sample_id) {
            out <- calculateUnivariate(x, type, features, colGraphName, s,
                exprs_values, BPPARAM, zero.policy,
                returnDF = TRUE,
                include_self = include_self, p.adjust.method = p.adjust.method,
                ...
            )
            local <- .is_local(type)
            if (local) {
                if (length(features) == 1L) {
                    names(out) <- features
                }
                localResults(x, type, features, sample_id = s) <- out
            } else {
                out <- .add_name_sample_id(out, s)
                features <- .symbol2id(x, features, swap_rownames)
                rowData(x)[features, names(out)] <- out
            }
        }
        x
    }
    if (is.null(type)) {
        fun_use
    } else {
        function(x, features = NULL, colGraphName = 1L, sample_id = "all",
                 exprs_values = "logcounts", BPPARAM = SerialParam(),
                 swap_rownames = NULL,
                 zero.policy = NULL, include_self = FALSE,
                 p.adjust.method = "BH", ...) {
            fun_use(
                x, type, features, colGraphName, sample_id,
                exprs_values, BPPARAM, swap_rownames, zero.policy, include_self,
                p.adjust.method, ...
            )
        }
    }
}
