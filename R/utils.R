.check_features <- function(x, features, colGeometryName = NULL) {
  # Check if features are in the gene count matrix or colData.
  # If not found, then assume that they're in the colGeometry
  if (is.null(features)) features <- rownames(x)
  features_assay <- intersect(features, rownames(x))
  features_coldata <- intersect(features, names(colData(x)))
  if (is.null(colGeometryName)) {
    features_colgeom <- NULL
  } else {
    cg <- colGeometry(x, type = colGeometryName)
    features_colgeom <- intersect(features, names(st_drop_geometry(cg)))
  }
  out <- list(assay = features_assay,
              coldata = features_coldata,
              colgeom = features_colgeom)
  if (all(lengths(out) == 0L)) {
    stop("None of the features are found in the SFE object.")
  }
  return(out)
}

.drop_null_list <- function(l) {
  null_inds <- vapply(l, is.null, FUN.VALUE = logical(1L))
  l[!null_inds]
}

.check_sample_id <- SpatialFeatureExperiment:::.check_sample_id

.is_discrete <- function(m) is.character(m) | is.factor(m) | is.logical(m)

.add_name_sample_id <- function(x, out, sample_id) {
  if (length(sampleIDs(x)) > 1L) {
    names(out) <- paste(names(out), sample_id, sep = "_")
  }
  out
}

#' @importFrom S4Vectors make_zero_col_DFrame
.initialize_featureData <- function(df) {
  if (is.null(attr(df, "featureData"))) {
    fd <- make_zero_col_DFrame(nrow = ncol(df))
    rownames(fd) <- colnames(df)
    attr(df, "featureData") <- fd
  }
  df
}

.add_fd <- function(x, df, res, features, sample_id, to_df_fun, name, to_df_params) {
  args_use <- c(list(out = res, name = name), to_df_params)
  res <- do.call(to_df_fun, args_use)
  res <- .add_name_sample_id(x, res, sample_id)
  df <- .initialize_featureData(df)
  fd <- attr(df, "featureData")
  fd[features, names(res)] <- res
  attr(df, "featureData") <- fd
  df
}

.get_feature_values <- function(sfe, features, sample_id,
                                colGeometryName = NULL) {
  features_list <- .check_features(sfe, features, colGeometryName)
  values <- list()
  sample_id_ind <- colData(sfe)$sample_id %in% sample_id
  if (!is.null(features_list[["assay"]])) {
    values_assay <- assay(sfe, exprs_values)[features_list[["assay"]],
                                             sample_id_ind, drop = FALSE]
    values_assay <- as.data.frame(as.matrix(t(values_assay)))
    values[["assay"]] <- values_assay
  }
  if (!is.null(features_list[["coldata"]]))
    values[["coldata"]] <- as.data.frame(colData(sfe)[sample_id_ind,
                                                      features_list[["coldata"]],
                                                      drop = FALSE])
  if (!is.null(features_list[["colgeom"]])) {
    cg <- colGeometry(sfe, colGeometryName, sample_id)
    values[["colgeom"]] <- st_drop_geometry(cg)[sample_id_ind,
                                                features_list[["colgeom"]],
                                                drop = FALSE]
  }
  if (length(values) > 1L) values <- do.call(cbind, values)
  values
}
