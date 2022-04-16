# Toy data to unit test functions in moran-geary.R
library(SpatialFeatureExperiment)
library(tidyverse)
library(Matrix)
data("visium_row_col")
# First sample
coords1 <- visium_row_col %>%
  filter(col < 6, row < 6)
coords1$row <- coords1$row * sqrt(3)

set.seed(29)
col_inds <- sample(1:13, 13)
row_inds <- sample(1:5, 13, replace = TRUE)
values <- sample(1:5, 13, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
colnames(mat) <- coords1$barcode

# Second sample
coords2 <- visium_row_col %>%
  filter(between(col, 8, 13), between(row, 8, 12))
coords2$row <- coords2$row * sqrt(3)

col_inds <- sample(1:15, 15)
row_inds <- sample(1:5, 15, replace = TRUE)
values <- sample(1:5, 15, replace = TRUE)
mat2 <- sparseMatrix(i = row_inds, j = col_inds, x = values)
colnames(mat2) <- coords2$barcode

rownames(mat) <- rownames(mat2) <- sample(LETTERS, 5)

sfe1 <- SpatialFeatureExperiment(list(counts = mat), colData = coords1,
                                 spatialCoordsNames = c("col", "row"),
                                 spotDiameter = 0.7)
sfe2 <- SpatialFeatureExperiment(list(counts = mat2), colData = coords2,
                                 spatialCoordsNames = c("col", "row"),
                                 sample_id = "sample02",
                                 spotDiameter = 0.7)
colGraph(sfe1, "visium1") <- findVisiumGraph(sfe1)
colGraph(sfe2, "visium2") <- findVisiumGraph(sfe2)
sfe <- cbind(sfe1, sfe2)
saveRDS(sfe, "inst/testdata/sfe.rds")
