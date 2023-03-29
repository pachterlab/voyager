library(SFEData)
library(spdep)
library(adespatial)
library(sf)

sfe <- McKellarMuscleData("small")
g <- findVisiumGraph(sfe)

test_that("listw2sparse gives correct results", {
    mat <- listw2sparse(g)
    expect_s4_class(mat, "dgCMatrix")
    expect_equal(nrow(mat), ncol(sfe))
    expect_equal(ncol(mat), ncol(sfe))
    expect_equal(Matrix::rowSums(mat > 0), card(g$neighbours), ignore_attr = TRUE)
    m2 <- listw2mat(g)
    expect_equal(as.matrix(mat), m2, ignore_attr = TRUE)
    expect_equal(rownames(mat), rownames(m2))
    expect_equal(rownames(mat), colnames(mat))
})

# Add a singleton to g
g_single <- g
g_single$neighbours <- c(g_single$neighbours, 0L)
class(g_single$neighbours) <- "nb"
attr(g_single, "region.id") <- c(attr(g_single, "region.id"), "foo")
g_single$weights <- c(g_single$weights, list(NULL))

test_that("Deal with singletons in listw2sparse", {
    mat <- listw2mat(g_single)
    n <- length(g_single$neighbours)
    expect_equal(nrow(mat), n)
    expect_equal(ncol(mat), n)
    expect_equal(Matrix::rowSums(mat)[n], 0, ignore_attr = TRUE)
})

test_that("moranBounds gives correct results", {
    # W scheme, default
    expect_equal(moranBounds(g), moran.bounds(g))
    # other schemes
    g2 <- findVisiumGraph(sfe, style = "B")
    expect_equal(moranBounds(g2), moran.bounds(g2))
    g3 <- findVisiumGraph(sfe, style = "S")
    expect_equal(moranBounds(g3), moran.bounds(g3))
})

nb1 <- grid2nb(d = c(5,5))
nb2 <- grid2nb(d = c(3,3))
attr(nb1, "region.id") <- LETTERS[1:25]
attr(nb2, "region.id") <- letters[1:9]
l1 <- nb2listw(nb1)
l2 <- nb2listw(nb2)
listws <- list(l1, l2)
names_expect <- c(LETTERS[1:25], letters[1:9])
test_that("Convert list of listws to one adjacency matrix", {
    mat <- multi_listw2sparse(listws)
    expect_s4_class(mat, "dgCMatrix")
    l_expect <- length(nb1) + length(nb2)
    expect_equal(nrow(mat), l_expect)
    expect_equal(ncol(mat), l_expect)
    expect_equal(rownames(mat), names_expect)
    expect_equal(colnames(mat), names_expect)
    expect_equal(as.matrix(mat[1:25,1:25]), listw2mat(l1), ignore_attr = TRUE)
    expect_equal(as.matrix(mat[26:34,26:34]), listw2mat(l2), ignore_attr = TRUE)
    expect_equal(sum(mat[26:34, 1:25]), 0)
    expect_equal(sum(mat[1:25, 26:34]), 0)
})
