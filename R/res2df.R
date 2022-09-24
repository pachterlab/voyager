# Convert different types of results into data frames
# Input: list of results, each element is for each feature
# Output: For global results, a DataFrame each row of which is for a feature,
# to be added to rowData or featureData.
# For local results, the results for each feature is reformatted as a vector or
# a matrix or data frame if necessary,
# and a DataFrame or data.frame (for geometries) each column of which is for
# a feature is returned.
.res2df <- function(out, type, local = FALSE, use_geometry = FALSE, ...) {
    if (local) {
        fun_use <- if (type %in% c("localmoran", "localmoran_perm"))
            .localmoran2df
        else if (type %in% c("localG", "localG_perm")) .localG2df
        else if (type == "localC_perm") .localCperm2df
        else identity
        out <- fun_use(out)
        out <- .value2df(out, use_geometry)
    } else {
        fun_use <- if (type %in% c("moran", "geary")) .moran2df
        else if (type %in% c("moran.mc", "geary.mc", "sp.mantel.mc", "EBImoran.mc",
                             "lee.mc")) .mcsim2df
        else if (type %in% c("moran.test", "geary.test", "globalG.test"))
            .htest2df
        else if (type == "sp.correlogram") .correlogram2df
        else .other2df
        out <- fun_use(out, type, ...)
    }
    out
}

.moran2df <- function(out, name, ...) {
    out <- lapply(out, unlist, use.names = TRUE)
    out <- Reduce(rbind, res)
    if (!is.matrix(out)) out <- t(as.matrix(out))
    rownames(out) <- rownames(x)
    out <- DataFrame(out)
    names(out)[1] <- name
    out
}

.mcsim2df <- function(out, name, ...) {
    out <- lapply(out, function(o) {
        o$res <- I(list(o$res))
        DataFrame(unclass(o))
    })
    rns <- names(out)
    out <- Reduce(rbind, out)
    rownames(out) <- rns
    colnames(out) <- paste(name, colnames(out), sep = "_")
    out
}

.htest2df <- function(out, name, ...) {
    # When it's not also mc.sim
    names_use <- c("statistic", "p.value", "alternative", "data.name", "method")
    out <- lapply(out, function(o) {
        o <- c(o[names_use], as.list(o[["estimate"]]))
        DataFrame(unclass(o))
    })
    rns <- names(out)
    out <- Reduce(rbind, out)
    rownames(out) <- rns
    colnames(out) <- paste(name, colnames(out), sep = "_")
    out
}

.correlogram2df <- function(out, name, ...) {
    other_args <- list(...)
    if ("method" %in% names(other_args)) method <- other_args[["method"]]
    else method <- "corr" # spdep's default
    out <- lapply(out, function(o) {
        colnames(o$res) <- c(method, "expectation", "variance")
        o
    })
    out <- lapply(out, function(o) o$res)
    out_df <- DataFrame(res = I(out))
    names(out_df) <- name
    out_df
}

.other2df <- function(out, name, ...) {
    if (!is.atomic(out)) out <- I(out)
    out_df <- DataFrame(res = out)
    names(out_df) <- name
    out_df
}

.localmoran2df <- function(out) {
    features <- names(out)
    out <- lapply(out, function(o) {
        o1 <- as.data.frame(o)
        quadr <- attr(o, "quadr")
        I(cbind(o1, quadr))
    })
    out
}

.attrmat2df <- function(out, attr_name, type) {
    if (attr_name %in% names(attributes(out[[1]]))) {
        out <- lapply(out, function(o) {
            attr_mat <- attr(o, attr_name)
            attr_mat <- cbind(o, attr_mat)
            colnames(attr_mat)[1] <- type
        })
    } else out
}

.localG2df <- function(out) .attrmat2df(out, "internals", "localG")

.localCperm2df <- function(out) .attrmat2df(out, "pseudo-p", "localC")
