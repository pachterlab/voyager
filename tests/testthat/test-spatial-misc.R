library(SFEData)
library(spdep)

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
