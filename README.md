
<img src="https://github.com/pachterlab/Voyager/raw/documentation/vignettes/voyager.jpg" width="1024"/>

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/Voyager.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Voyager)
<!-- badges: end -->

> "Everything is related to everything else. But near things are more related than distant things." - Tobler's first law of geography 

[`Voyager`](https://bioconductor.org/packages/devel/bioc/html/Voyager.html) is a package that facilitates exploratory spatial data analysis and visualization for spatial genomics data represented by [`SpatialFeatureExperiment`](https://bioconductor.org/packages/devel/bioc/html/SpatialFeatureExperiment.html) objects. 

`Voyager` and `SpatialFeatureExperiment` were developed within the Bioconductor ecosystem, where they build on several existing objects and tools. Single cell RNA-seq data and metadata can be represented with [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) S4 class objects, and these can be utilized for exploratory data analysis and visualization using the [`scater`](https://bioconductor.org/packages/release/bioc/html/scater.html), [`scran`](https://bioconductor.org/packages/release/bioc/html/scran.html), or [`scuttle`](https://bioconductor.org/packages/release/bioc/html/scuttle.html) packages. The  
[`SpatialExperiment`](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) class extends `SingleCellExperiments` to allow for representation of spatial genomics data. `SpatialFeatureExperiment` extends `SpatialExperiment` with Simple Features from [`sf`](https://r-spatial.github.io/sf/). 

Thus, `Voyager` : `SpatialFeatureExperiment` is as `scater` / `scran` / `scuttle` : `SingleCellExperiment`. `Voyager` also builds on the geospatial tradition, especially the [`spdep`](https://r-spatial.github.io/spdep/) package, which is one of the main R packages for spatial dependence analyses. Specifically, `Voyager` focuses on spatial autocorrelation, which measures the extent of similarity or dissimilarity of spatially proximal regions, and that can be quantified in terms of length scale, and variation in space.

Questions, suggestions, or problems should be submitted as GitHub issues](https://github.com/pachterlab/Voyager/issues).

<img src="https://github.com/pachterlab/Voyager/raw/documentation/vignettes/voyager_schematics.png" width="1024"/>

<!--- About the banner: USS Voyager resting on N San Gabriel Canyon Rd, along north fork San Gabriel River, north of Glendora, LA county --->
