library(SFEData)
library(scater)

sfe <- McKellarMuscleData("small")
sfe <- sfe[,sfe$in_tissue]
sfe <- logNormCounts(sfe)
inds <- order(Matrix::rowSums(logcounts(sfe)), decreasing = TRUE)[1:50]
mat <- logcounts(sfe)[inds,]
g <- colGraph(sfe, "visium") <- findVisiumGraph(sfe)
test_that("Correct output structure of multispati_rsp", {
    out <- multispati_rsp(t(mat), listw = g, nfposi = 10, nfnega = 10)
    expect_equal(colnames(out), paste0("PC", 1:20))
    expect_equal(rownames(out), colnames(sfe))
    expect_true(is.numeric(out))
    loadings <- attr(out, "rotation")
    expect_equal(colnames(loadings), paste0("PC", 1:20))
    expect_equal(rownames(loadings), rownames(mat))
    expect_true(is.numeric(loadings))
    eigs <- attr(out, "eig")
    expect_equal(length(eigs), 20)
    expect_true(is.numeric(eigs))
    expect_true(all(eigs[1:10] > 0))
    expect_true(all(diff(eigs) < 0))
    expect_true(all(eigs[11:20] < 0))
})

ref <- multispati_rsp(t(mat), listw = g, nfposi = 10, nfnega = 10)
test_that("CalculateMultivariate for matrix", {
    out <- calculateMultivariate(mat, "multispati", listw = g,
                                 nfposi = 10, nfnega = 10)
    expect_equal(out, ref)
})

test_that("CalculateMultivariate for SFE, one sample", {
    out <- calculateMultivariate(sfe, "multispati", colGraphName = "visium",
                                 nfposi = 10, nfnega = 10, subset_row = inds)
    expect_equal(out, ref)
})

sfe2 <- McKellarMuscleData("small2")
sfe2 <- sfe2[,sfe2$in_tissue]
sfe2 <- logNormCounts(sfe2)
colGraph(sfe2, "visium") <- findVisiumGraph(sfe2)
sfe3 <- SpatialFeatureExperiment::cbind(sfe, sfe2)
# More general case when different samples don't occupy distinct blocks
set.seed(29)
rand_inds <- sample(seq_len(ncol(sfe3)), ncol(sfe3))
sfe3 <- sfe3[,rand_inds]
test_that("When results from different samples are concatenated", {
    out <- calculateMultivariate(sfe3, "multispati", colGraphName = "visium",
                                 sample_action = "separate",
                                 nfposi = 10, nfnega = 10, subset_row = inds)
    expect_equal(rownames(out), colnames(sfe3))
    expect_equal(colnames(out), paste0("PC", 1:20))
    expect_true(is.numeric(out))
    loadings <- attr(out, "rotation")
    expect_true(is.numeric(loadings))
    expect_equal(rownames(loadings), rownames(sfe3)[inds])
    expect_equal(colnames(loadings), paste0("PC", 1:20))
    expect_equal(dimnames(loadings)[[3]], sampleIDs(sfe3))
    eigs <- attr(out, "eig")
    expect_true(is.matrix(eigs))
    expect_equal(colnames(eigs), sampleIDs(sfe3))
    expect_equal(nrow(eigs), 20)
    expect_true(all(eigs[1:10,] > 0))
    expect_true(all(eigs[11:20,] < 0))
})

test_that("Correct output structure for joint", {
    out <- calculateMultivariate(sfe3, "multispati", colGraphName = "visium",
                                 sample_action = "joint",
                                 nfposi = 10, nfnega = 10, subset_row = inds)
    expect_equal(rownames(out), colnames(sfe3))
    expect_equal(colnames(out), paste0("PC", 1:20))
    expect_true(is.numeric(out))
    loadings <- attr(out, "rotation")
    expect_true(is.numeric(loadings))
    expect_true(is.matrix(loadings))
    expect_equal(colnames(loadings), paste0("PC", 1:20))
    expect_equal(rownames(loadings), rownames(sfe3)[inds])
    eigs <- attr(out, "eig")
    expect_equal(length(eigs), 20)
    expect_true(is.numeric(eigs))
    expect_true(all(eigs[1:10] > 0))
    expect_true(all(diff(eigs) < 0))
    expect_true(all(eigs[11:20] < 0))
})

test_that("Correctly add results to the SFE object", {
    sfe <- runMultivariate(sfe, "multispati", subset_row = inds, nfposi = 10,
                           nfnega = 10)
    out <- reducedDim(sfe, "multispati")
    expect_equal(out, ref)
})
# 6. Correct output structure when output is a vector
# including the message on reducedDim
