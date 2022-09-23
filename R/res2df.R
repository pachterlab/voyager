# Convert different types of results into data frames
# Input: list of results

.res2df <- function(out, type, ...) {
    fun_use <- if (type %in% c("moran", "geary")) .moran2df
    else if (type %in% c("moran.mc", "geary.mc", "sp.mantel.mc", "EBImoran.mc",
                         "lee.mc")) .mcsim2df
    else if (type %in% c("moran.test", "geary.test", "globalG.test"))
        .htest2df
    else if (type == "sp.correlogram") .correlogram2df
    else .other2df
    fun_use(out, type, ...)
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
