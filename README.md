
# From geospatial to spatial transcriptomics

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/Voyager.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Voyager)
[![codecov](https://codecov.io/github/pachterlab/Voyager/branch/master/graph/badge.svg?token=RCIXA7AQER)](https://codecov.io/github/pachterlab/Voyager)
<!-- badges: end -->

This package brings the tradition of geospatial statistics to spatial omics by wrapping classical geospatial packages such as `spdep` and `adespatial` to be used with the SpatialFeatureExperiment class, which extends SpatialExperiment with sf.

The companion [website](https://pachterlab.github.io/voyager/) for this package includes vignettes that showcase the functionality of `Voyager` in the context of Visium, Slide-seq V2, CosMx, Xenium, and MERFISH.  
## Installation

This package is in Bioconductor version 3.16, which is currently the devel version. Install with

```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("Voyager", version = "devel")
```

The development version of Voyager from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("pachterlab/Voyager")
```
