library(SFEData)
library(spdep)
library(adespatial)

sfe <- McKellarMuscleData("small")
g <- findVisiumGraph(sfe)

test_that("listw2sparse gives correct results", {
    mat <- listw2sparse(g)
    expect_s4_class(mat, "dgCMatrix")
    expect_equal(nrow(mat), ncol(sfe))
    expect_equal(ncol(mat), ncol(sfe))
    expect_equal(Matrix::rowSums(mat > 0), card(g$neighbours))
    m2 <- listw2mat(g)
    dimnames(m2) <- NULL
    expect_equal(as.matrix(mat), m2)
})

test_that("moranBounds gives correct results", {
    # W scheme, default
    expect_equal(moranBounds(g), moran.bounds(g))
    # other schemes
    g2 <- findVisiumGraph(sfe, style = "B")
    expect_equal(moranBounds(g2), moran.bounds(g2))
    g3 <- findVisiumGraph(sfe, style = "C")
    expect_equal(moranBounds(g3), moran.bounds(g3))
})
