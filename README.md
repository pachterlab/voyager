
# From geospatial to spatial transcriptomics

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/Voyager.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Voyager)
[![codecov](https://codecov.io/gh/pachterlab/voyager/branch/main/graph/badge.svg?token=RCIXA7AQER)](https://codecov.io/gh/pachterlab/voyager)
<!-- badges: end -->

This package brings the tradition of geospatial statistics to spatial omics by wrapping classical geospatial packages such as `spdep` and `adespatial` to be used with the SpatialFeatureExperiment class, which extends SpatialExperiment with sf.

The [companion website](https://pachterlab.github.io/voyager/) for this package includes vignettes that showcase the functionality of `Voyager` in the context of the Visium, Slide-seq V2, CosMx, Xenium, and MERFISH technologies.  
## Installation

This package is in Bioconductor version 3.16 and above. Install with

```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(version = "3.16") # Or a higher version in the future
BiocManager::install("Voyager")
```

The main branch in this repo is the release version. The development version of Voyager can be from [GitHub](https://github.com/) with:

```r
# install.packages("remotes")
remotes::install_github("pachterlab/voyager", ref = "devel")
```

Or from Bioconductor with:

```r
BiocManager::install("Voyager", version = "devel")
```
