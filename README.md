
<img src="https://github.com/pachterlab/Voyager/raw/documentation/vignettes/voyager.jpg" width="1024"/>
 P
<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/Voyager.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Voyager)
<!-- badges: end -->

<img src="https://github.com/pachterlab/Voyager/raw/documentation/vignettes/voyager_schematics.png" width="1024"/>


> "Everything is related to everything else. But near things are more related than distant things." - Tobler's first law of geography 

Within the Bioconductor ecosystem, single cell RNA-seq data and metadata can be represented with [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) S4 class objects for exploratory data analysis and visualization using the [`scater`](https://bioconductor.org/packages/release/bioc/html/scater.html), [`scran`](https://bioconductor.org/packages/release/bioc/html/scran.html), or [`scuttle`](https://bioconductor.org/packages/release/bioc/html/scuttle.html) packages. The  
[`SpatialExperiment`](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) class extends `SingleCellExperiments` to allow for representation of spatial genomics data. 

We have developed the [`SpatialFeatureExperiment`](https://bioconductor.org/packages/devel/bioc/html/SpatialFeatureExperiment.html) S4 class, which extends `SpatialExperiment`, to provide infrastructure for vector spatial data, and the companion [`Voyager`](https://bioconductor.org/packages/devel/bioc/html/Voyager.html) package facilitates exploratory data analysis and visualization for `SpatialFeatureExperiment` objects. Thus, `Voyager` : `SpatialFeatureExperiment` is as like `scater` / `scran` / `scuttle` : `SingleCellExperiment`. `Voyager` also builds on the geospatial tradition, especially the [`spdep`](https://r-spatial.github.io/spdep/) package, which is one of the main R packages for spatial dependence analyses. Specifically, `Voyager` focuses on spatial autocorrelation, which measures the extent of similarity or dissimilarity of spatially proximal regions, and that can be quantified in terms of length scale, and variation in space.

Please open a [GitHub issue](https://github.com/pachterlab/Voyager/issues) if you have questions, suggestions, or have encountered problems.

<!--- About the banner: USS Voyager resting on N San Gabriel Canyon Rd, along north fork San Gabriel River, north of Glendora, LA county --->
