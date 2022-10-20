
# From geospatial to spatial transcriptomics

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Rules of the `documentation` branch:

* The purpose of this branch is to make a `pkgdown` website with longer and detailed vignettes that would make the installed size of this package way too large to comply with Bioconductor's 5 MB rule.
* Don't change anything outside the `vignettes` directory in this branch. If the code doesn't work, change in the `main` or `devel` branch and merge into this branch. This way the large vignettes won't get into the `main` branch and the code is kept consistent, which is important since the `pkgdown` website also documents all the functions of this package.
* Exception to the previous rule: you may add packages to the Suggests field in DESCRIPTION for extra packages used in the vignettes.
* To trigger a new GitHub Action build of the website, add a tag with `git tag <tag name>` and push it. The Action will only be run if a tag is pushed. The tag name can be something like "update_visium".
* The file `vignettes/ref.bib` is automatically synced from Paperpile. Don't edit by hand.

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
