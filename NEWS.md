# Version 1.1.10 (03/07/2023)
* Record parameters used to get spatial results
* Force users to use a new name when running the same method with different parameters

# Version 1.1.9 (02/12/2023)
* Deprecated show_symbol argument, replacing with swap_rownames to be consistent with scater

# Version 1.1.7
* Added bbox argument to spatial plotting functions to zoom in with a bounding box

# Version 1.0.10 (02/23/2023)
* Added plotColDataFreqpoly when the y axis needs to be log transformed. That doesn't work with stacked histograms and using position = "identity" causes some bars to be covered.

# Version 1.0.9
* Fixed the bug of hardcoded ncol in plotDimLoadings.

# Version 1.0.8
* Flipped the divergent palettes so the warm color means high value.

# Version 1.0.7
* Fixed bug with assigning local results for each sample for colData, colGeometry, and annotGeometry.

# Version 1.0.5
* Removed aes_string(), which is deprecated.
* Fixed bug when show_symbol = TRUE and "symbol" column is absent from rowData.

# Version 1.0.0

* First version on Bioconductor
* Univariate local and global spatial statistics based on spdep
* Plotting functions: gene expression and metadata in space, results of local spatial analyses, plot dimension reductions in space, and plot correlograms and Monte Carlo simulation results
