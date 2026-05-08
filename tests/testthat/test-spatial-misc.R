library(SFEData)
library(spdep)
library(sf)

sfe <- McKellarMuscleData("small")
g <- findVisiumGraph(sfe)

# Copied and pasted code from adespatial, which got kicked out of CRAN
bicenter.wt <- function (X, row.wt = rep(1, nrow(X)), col.wt = rep(1, ncol(X))) {
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    if (length(row.wt) != n)
        stop("length of row.wt must equal the number of rows in x")
    if (any(row.wt < 0) || (sr <- sum(row.wt)) == 0)
        stop("weights must be non-negative and not all zero")
    row.wt <- row.wt/sr
    if (length(col.wt) != p)
        stop("length of col.wt must equal the number of columns in x")
    if (any(col.wt < 0) || (st <- sum(col.wt)) == 0)
        stop("weights must be non-negative and not all zero")
    col.wt <- col.wt/st
    row.mean <- apply(row.wt * X, 2, sum)
    col.mean <- apply(col.wt * t(X), 2, sum)
    col.mean <- col.mean - sum(row.mean * col.wt)
    X <- sweep(X, 2, row.mean)
    X <- t(sweep(t(X), 2, col.mean))
    return(X)
}

moran.bounds <- function(listw){
    W <- listw2mat(listw)
    W <- (W + t(W)) / 2
    res <- nrow(W) / sum(W) *  range(eigen(bicenter.wt(W), only.values = TRUE))
    names(res) <- c("Imin", "Imax")
    return(res)
}

test_that("moranBounds gives correct results", {
    # W scheme, default
    expect_equal(moranBounds(g), moran.bounds(g))
    # other schemes
    g2 <- findVisiumGraph(sfe, style = "B")
    expect_equal(moranBounds(g2), moran.bounds(g2))
    g3 <- findVisiumGraph(sfe, style = "S")
    expect_equal(moranBounds(g3), moran.bounds(g3))
})
