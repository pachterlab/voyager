library(SFEData)
library(scater)

sfe <- McKellarMuscleData("small")
sfe <- sfe[,sfe$in_tissue]
sfe <- logNormCounts(sfe)
inds <- order(Matrix::rowSums(logcounts(sfe)), decreasing = TRUE)[1:50]
mat <- logcounts(sfe)[inds,]
g <- colGraph(sfe, "visium") <- findVisiumGraph(sfe)

test_that("Error when a univariate method is used", {
    expect_error(calculateMultivariate(mat, type = "moran", listw = g),
                 "`type` must be a multivariate method.")
})

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
    expect_message({
        sfe <- runMultivariate(sfe, "multispati", subset_row = inds, nfposi = 10,
                               nfnega = 10, dest = "colData")
    }, "Matrix or array outputs can only be stored in reducedDims.")
})

test_that("When output is a vector", {
    out <- calculateMultivariate(sfe, "localC_multi", colGraphName = "visium",
                                 subset_row = inds)
    expect_vector(out)
    sfe <- runMultivariate(sfe, "localC_multi", colGraphName = "visium",
                           subset_row = inds, dest = "colData")
    expect_true("localC_multi" %in% names(colData(sfe)))
    expect_equal(sfe$localC_multi, out)
    expect_message({
        sfe <- runMultivariate(sfe, "localC_multi", colGraphName = "visium",
                               subset_row = inds, dest = "reducedDim")
    }, "Vector output can only be stored in colData.")
})

names_expect_lc <- c(
    "localC_perm_multi", "E.Ci", "Var.Ci", "Z.Ci", "Pr(z != E(Ci))",
    "Pr(z != E(Ci)) Sim", "Pr(folded) Sim", "Skewness",
    "Kurtosis", "-log10p Sim", "-log10p_adj Sim", "cluster"
)
test_that("When output is a data frame", {
    df <- calculateMultivariate(sfe, "localC_perm_multi", colGraphName = "visium",
                                subset_row = inds, nsim = 99)
    expect_s3_class(df, "data.frame")
    expect_equal(names(df), names_expect_lc)
    expect_equal(nrow(df), ncol(sfe))
})

test_that("Add data frame output to colData", {
    sfe <- runMultivariate(sfe, "localC_perm_multi", colGraphName = "visium",
                           subset_row = inds, nsim = 99, dest = "colData")
    names_expect <- names_expect_lc
    names_expect[-1] <- paste("localC_perm_multi", names_expect[-1], sep = "_")
    expect_true(all(names_expect %in% names(colData(sfe))))
})

test_that("Add data frame output to reducedDim", {
    sfe <- runMultivariate(sfe, "localC_perm_multi", colGraphName = "visium",
                           subset_row = inds, nsim = 99, dest = "reducedDim")
    df <- reducedDim(sfe, "localC_perm_multi")
    expect_s3_class(df, "data.frame")
    expect_equal(names(df), names_expect_lc)
})
