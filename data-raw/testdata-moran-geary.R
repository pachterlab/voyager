# Toy data to unit test functions in moran-geary.R
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(sf)
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
colGraph(sfe1, "visium") <- findVisiumGraph(sfe1)
colGraph(sfe2, "visium") <- findVisiumGraph(sfe2)
sfe <- cbind(sfe1, sfe2)

# Add annotGeometry and annotGraph
set.seed(29)
bbox1 <- st_bbox(spotPoly(sfe, "sample01"))
annot_x <- runif(5, min = bbox1["xmin"], max = bbox1["xmax"])
annot_y <- runif(5, min = bbox1["ymin"], max = bbox1["ymax"])
annot_coords <- data.frame(x = annot_x, y = annot_y, sample_id = "sample01")
annot_geom <- df2sf(annot_coords)
annot_geom <- st_buffer(annot_geom, runif(5, max = 0.5))
annotGeometry(sfe, "annot", "sample01") <- annot_geom
annotGraph(sfe, "annot_tri", "sample01") <-
  findSpatialNeighbors(sfe, "sample01", type = "annot", MARGIN = 3)

# Add geometry metadata
set.seed(29)
colGeometry(sfe, "spotPoly", "all")$foo <- rnorm(ncol(sfe))
colData(sfe)$nCounts <- colSums(counts(sfe))
annotGeometry(sfe, "annot", "all")$bar <- rnorm(5)
colGeometry(sfe, "centroids", "all") <- suppressWarnings(st_centroid(spotPoly(sfe, "all")))
colGeometry(sfe, "centroids", "all")$category <- sample(LETTERS[22:26], ncol(sfe), replace = TRUE)

saveRDS(sfe, "inst/testdata/sfe.rds")
