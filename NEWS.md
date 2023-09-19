# Version 1.3.1 (05/15/2023)
* Removed functions and arguments deprecated in 1.2.0

# Version 1.2.5 (08/18/2023)
* Use imgRaster getter rather than the S4 no-no of @image to get images to plot, as the latter will no longer work as of SFE 1.2.3 that wraps SpatRaster images when saving RDS. Reading RDS won't unwrap so images need to be unwrapped when they're needed.

>>>>>>> main
# Version 1.2.4 (07/04/2023)
* Remove useNames = NA warning when calling MULTISPATI; the warning comes from
generic of colVars.
* Use algebraic eigenvalues for MULTISPATI when either nfposi or nfnega is 0
* Added bins_contour argument to moranPlot to change the number of bins in cell
density contours

# Version 1.2.3 (05/04/2023)
* Fix bug when plotting a feature with illegal name alongside another feature
with legal name
* Make sure runBivariate and calculateBivariate use gene symbols in results even
if Ensembl IDs are specified when swap_rownames is set
* Change secondary sequential palette in the light theme to YlOrRd so it's more
distinguishable from the Blues primary palette at low values

# Version 1.2.2 (04/26/2023)
* Some minor bugs: runBivariate gets correct feature names when only feature1 is
specified and swap_rownames is used to show gene symbol
* Correct output for cross variogram maps for only one pair of genes
* Added default_attr to localmoran_bv's SFEMethod
* Don't plot attribute when localResult is a vector and there's no default attr
* When plotting multiple features, the panels follow the same order the features
are specified
* Allow illegal characters in names of colData and reducedDims in plots
* Plot only one component in spatialReducedDim with the components argument
* Deprecate plotColDataBin2D and plotRowDataBin2D

# Version 1.1.12 (04/22/2023)
* Plot image behind geometries in all functions that plot geometries
* Added dark theme support for functions that plot geometries

# Version 1.1.11 (04/05/2023)
* Added MULTISPATI PCA
* Added multivariate local Geary's C from Anselin 2019
* Added calculateMultivariate as a unified user interface to multivariate spatial analyses
* Variogram and variogram map with gstat and related plotting functions
* Allow non-standard names for local results in plotLocalResult

# Version 1.1.10 (03/07/2023)
* Record parameters used to get spatial results
* Force users to use a new name when running the same method with different parameters

# Version 1.1.9 (02/12/2023)
* Deprecated show_symbol argument, replacing with swap_rownames to be consistent with scater

# Version 1.1.7
* Added bbox argument to spatial plotting functions to zoom in with a bounding box

# Version 1.0.10 (02/23/2023)
* Added plotColDataFreqpoly when the y axis needs to be log transformed. That doesn't work with stacked histograms and using position = "identity" causes some bars to be covered.

# Version 1.0.9 (02/03/2023)
* Fixed the bug of hardcoded ncol in plotDimLoadings.

# Version 1.0.8 (01/26/2023)
* Flipped the divergent palettes so the warm color means high value.

# Version 1.0.7 (01/11/2023)
* Fixed bug with assigning local results for each sample for colData, colGeometry, and annotGeometry.

# Version 1.0.5 (12/02/2022)
* Removed aes_string(), which is deprecated.
* Fixed bug when show_symbol = TRUE and "symbol" column is absent from rowData.

# Version 1.0.0 (11/02/2022)

* First version on Bioconductor
* Univariate local and global spatial statistics based on spdep
* Plotting functions: gene expression and metadata in space, results of local spatial analyses, plot dimension reductions in space, and plot correlograms and Monte Carlo simulation results
