.check_features_all <- function(x, features, colGeometryName) {
  # Check if features are in the gene count matrix or colData.
  # If not found, then assume that they're in the colGeometry
  features_assay <- intersect(features, rownames(x))
  features_coldata <- intersect(features, names(colData(x)))
  if (missing(colGeometryName)) {
    features_colgeom <- NULL
  } else {
    cg <- colGeometry(x, type = colGeometryName)
    features_colgeom <- intersect(features, names(st_drop_geometry(cg)))
  }
  out <- list(assay = featuers_assay,
              coldata = features_coldata,
              colgeom = features_colgeom)
  if (all(lengths(out) == 0L)) {
    stop("None of the features are found in the SFE object.")
  }
  return(out)
}

.drop_null_list <- function(l) as.list(unlist(l, use.names = TRUE))
