# Convert different types of results into data frames
# Input: list of results, each element is for each feature
# Output: For global results, a DataFrame each row of which is for a feature,
# to be added to rowData or featureData.
# For local results, the results for each feature is reformatted as a DataFrame,
# and a DataFrame or data.frame (for geometries) each column of which is for
# a feature is returned.
#' @importFrom spdep p.adjustSP
.add_log_p <- function(out, nb, p.adjust.method) {
  if (!is.vector(out)) {
    ind_sim <- grepl("^Pr\\(z.+Sim$", colnames(out))
    ind <- grepl("^Pr\\(", colnames(out))
    if (any(ind_sim)) {
      col_use <- colnames(out)[ind_sim][1]
      out <- cbind(out, `-log10p Sim` = -log10(out[,col_use]),
                   `-log10p_adj Sim` = -log10(p.adjustSP(out[,col_use],nb,
                                                        method = p.adjust.method)))
    } else if (any(ind)) {
      col_use <- colnames(out)[ind][1]
      out <- cbind(out, `-log10p` = -log10(out[,col_use]),
                   `-log10p_adj` = -log10(p.adjustSP(out[,col_use], nb,
                                                    method = p.adjust.method)))
    }
  }
  out
}

.moran2df <- function(out, name, ...) {
    rns <- names(out)
    out <- lapply(out, unlist, use.names = TRUE)
    out <- Reduce(rbind, out)
    if (!is.matrix(out)) out <- t(as.matrix(out))
    rownames(out) <- rns
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
        o[names_use] <- lapply(
            o[names_use],
            function(n) {
                if (is.matrix(n) && all(dim(n) == 1)) {
                    n <- as.vector(n)
                }
                n
            }
        )
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
    if ("method" %in% names(other_args)) {
        method <- other_args[["method"]]
    } else {
        method <- "I"
    }
    name <- paste(name, method, sep = "_")
    if (method %in% c("I", "C")) {
        out <- lapply(out, function(o) {
            colnames(o$res) <- c(method, "expectation", "variance")
            o
        })
    }
    out <- lapply(out, function(o) o$res)
    out_df <- DataFrame(res = I(out))
    names(out_df) <- name
    rownames(out_df) <- names(out)
    out_df
}

.other2df <- function(out, name, ...) {
    if (!is.atomic(out)) out <- I(out)
    out_df <- DataFrame(res = out)
    names(out_df) <- name
    rownames(out_df) <- names(out)
    out_df
}

.localmoran2df <- function(out, nb, p.adjust.method) {
    lapply(out, function(o) {
        o1 <- as.data.frame(o)
        quadr <- attr(o, "quadr")
        I(.add_log_p(cbind(o1, quadr), nb, p.adjust.method))
    })
}

.attrmat2df <- function(out, attr_name, type, nb, p.adjust.method) {
    if (attr_name %in% names(attributes(out[[1]]))) {
        lapply(out, function(o) {
            attr_mat <- attr(o, attr_name)
            attr_mat <- cbind(o, attr_mat)
            colnames(attr_mat)[1] <- type
            res <- .add_log_p(attr_mat, nb, p.adjust.method)
            if ("cluster" %in% names(attributes(o))) {
                res <- as.data.frame(res)
                res$cluster <- attr(o, "cluster")
            }
            res
        })
    } else {
        lapply(out, as.vector)
    }
}

.localG2df <- function(out, nb, p.adjust.method)
    .attrmat2df(out, "internals", "localG", nb, p.adjust.method)

.localCperm2df <- function(out, nb, p.adjust.method)
    .attrmat2df(out, "pseudo-p", "localC", nb, p.adjust.method)

.localCpermmulti2df <- function(out, nb, p.adjust.method) {
    .attrmat2df(list(out), "pseudo-p", "localC_perm_multi", nb, p.adjust.method)[[1]]
}

.LOSHmc2df <- function(out, nb, p.adjust.method)
    lapply(out, .add_log_p, nb = nb, p.adjust.method = p.adjust.method)

.to_df_identity <- function(out, nb, p.adjust.method) out
