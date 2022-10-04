
# From geospatial to spatial transcriptomics

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package brings the tradition of geospatial statistics to spatial omics by wrapping classical geospatial packages such as `spdep` and `adespatial` to be used with the SpatialFeatureExperiment class, which extends SpatialExperiment with sf.

## Installation

This package has been submitted to Bioconductor. Once accepted, it can be installed with 

```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("Voyager")
```

The development version of Voyager from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("pachterlab/Voyager")
```
