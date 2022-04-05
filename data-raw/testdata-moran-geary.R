# Toy data to unit test functions in moran-geary.R
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
sfe_visium <- readRDS(system.file("testdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
g_visium <- readRDS(system.file("testdata/colgraph_visium.rds",
                                package = "SpatialFeatureExperiment"))
g_visium2 <- readRDS(system.file("testdata/colgraph_visium2.rds",
                                 package = "SpatialFeatureExperiment"))
