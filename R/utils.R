#' @importFrom SpatialFeatureExperiment .value2df .check_features
#'   .warn_symbol_duplicate .symbol2id .check_sample_id .rm_empty_geometries

.drop_null_list <- function(l) {
  null_inds <- vapply(l, is.null, FUN.VALUE = logical(1L))
  l[!null_inds]
}

.is_discrete <- function(m) is.character(m) | is.factor(m) | is.logical(m)

.add_sample_id <- function(name, sample_id) {
  if (!grepl(sample_id, name)) {
    name <- paste(name, sample_id, sep = "_")
  }
  name
}

#' @importFrom SpatialFeatureExperiment sampleIDs annotGeometry annotGeometry<-
#' colGeometry colGeometry<-
#' @importFrom SummarizedExperiment colData<-
#' @importFrom methods is
.add_name_sample_id <- function(out, sample_id) {
  names(out) <- vapply(names(out), .add_sample_id, sample_id = sample_id,
                       FUN.VALUE = character(1))
  out
}

#' @importFrom S4Vectors make_zero_col_DFrame
#' @importFrom SingleCellExperiment int_metadata int_metadata<-
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
  res <- .add_name_sample_id(res, sample_id)
  df <- .initialize_featureData(df)
  fd <- attr(df, "featureData")
  fd[features, names(res)] <- res
  attr(df, "featureData") <- fd
  df
}

.initialize_fd_dimData <- function(x, MARGIN) {
  fd_name <- switch(MARGIN, "rowFeatureData", "colFeatureData")
  if (is.null(int_metadata(x)[[fd_name]])) {
    rownames_use <- switch(MARGIN, names(rowData(x)), names(colData(x)))
    fd <- make_zero_col_DFrame(nrow = length(rownames_use))
    rownames(fd) <- rownames_use
    int_metadata(x)[[fd_name]] <- fd
  }
  x
}

# Because adding a new column to S4 DataFrame will remove the attributes
# Put the featureData of colData and rowData in int_metadata instead
.add_fd_dimData <- function(x, MARGIN, res, features, sample_id, type, ...) {
  res <- .res2df(res, type, ...)
  res <- .add_name_sample_id(res, sample_id)
  x <- .initialize_fd_dimData(x, MARGIN)
  fd_name <- switch(MARGIN, "rowFeatureData", "colFeatureData")
  fd <- int_metadata(x)[[fd_name]]
  fd[features, names(res)] <- res
  int_metadata(x)[[fd_name]] <- fd
  x
}

#' Get metadata of colData and rowData
#'
#' Results of spatial analyses on columns in \code{colData} and \code{rowData}
#' are stored in \code{int_metadata(sfe)}, or internal metadata. This function
#' allows the users to access these results.
#'
#' @param sfe An SFE object.
#' @return A \code{DataFrame}.
#' @export
#' @name colFeatureData
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Moran's I for colData
#' sfe <- colDataMoransI(sfe, "nCounts")
#' colFeatureData(sfe)
colFeatureData <- function(sfe) {
  int_metadata(sfe)$colFeatureData
}

#' @rdname colFeatureData
#' @export
rowFeatureData <- function(sfe) {
  int_metadata(sfe)$rowFeatureData
}

.get_feature_values <- function(sfe, features, sample_id,
                                colGeometryName = NULL, annotGeometryName = NULL,
                                exprs_values = "logcounts", cbind_all = TRUE,
                                show_symbol = TRUE) {
    features_list <- .check_features(sfe, features, colGeometryName, annotGeometryName)
    values <- list()
    sample_id_ind <- colData(sfe)$sample_id %in% sample_id
    if (length(features_list[["assay"]])) {
        values_assay <- assay(sfe, exprs_values)[features_list[["assay"]],
                                                 sample_id_ind, drop = FALSE]
        # So symbol is shown instead of Ensembl ID
        if ("symbol" %in% names(rowData(sfe)) && show_symbol) {
            rownames(values_assay) <- rowData(sfe)[rownames(values_assay), "symbol"]
        }
        values_assay <- as.data.frame(as.matrix(t(values_assay)))
        values[["assay"]] <- values_assay
    }
    if (length(features_list[["coldata"]]))
        values[["coldata"]] <- as.data.frame(colData(sfe)[sample_id_ind,
                                                          features_list[["coldata"]],
                                                          drop = FALSE])
    if (length(features_list[["colgeom"]])) {
        cg <- colGeometry(sfe, colGeometryName, sample_id)
        values[["colgeom"]] <- st_drop_geometry(cg)[sample_id_ind,
                                                    features_list[["colgeom"]],
                                                    drop = FALSE]
    }
    if (length(features_list[["annotgeom"]])) {
        ag <- annotGeometry(sfe, annotGeometryName, sample_id)
        sample_id_ind2 <- ag$sample_id %in% sample_id
        values[["annotgeom"]] <- st_drop_geometry(ag)[sample_id_ind,
                                                      features_list[["annotgeom"]],
                                                      drop = FALSE]
    }
    if (cbind_all) {
        if (length(values) > 1L) {
            if (length(features_list[["annotgeom"]]))
                warning("annotGeometry values cannot be cbinded to other values.")
            values <- do.call(cbind, values[c("assay", "coldata", "colgeom")])
        } else values <- values[[1]]
    }
    values
}

.is_na_list <- function(l) {
  vapply(l, function(d) isTRUE(is.na(d)), logical(1))
}

.get_not_na_items <- function(df, features, colname_use) {
  if (is(df, "sf")) df <- st_drop_geometry(df)
  if (!colname_use %in% names(df)) return(NULL)
  out <- setNames(df[features, colname_use], features)
  out[!.is_na_list(out)]
}

.get_feature_metadata <- function(sfe, features, name, sample_id,
                                  colGeometryName, annotGeometryName,
                                  show_symbol) {
  colname_use <- .add_sample_id(name, sample_id)
  out_rd <- out_cd <- out_cg <- out_ag <- NULL
  features_rd <- intersect(features, rownames(sfe))
  if (!length(features_rd) && "symbol" %in% names(rowData(sfe))) {
    features_symbol <- intersect(features, rowData(sfe)$symbol)
    .warn_symbol_duplicate(sfe, features_rd)
    features <- setdiff(features, features_symbol)
    features_rd <- rownames(sfe)[match(features_symbol, rowData(sfe)$symbol)]
    if (all(is.na(features_rd))) features_rd <- NULL
  }
  if (length(features_rd)) {
    out_rd <- .get_not_na_items(rowData(sfe), features_rd, colname_use)
    features <- setdiff(features, names(out_rd))
    if ("symbol" %in% names(rowData(sfe)) && show_symbol) {
      names(out_rd) <- rowData(sfe)[names(out_rd), "symbol"]
      .warn_symbol_duplicate(sfe, names(out_rd))
    }
  }
  features_cd <- intersect(features, names(colData(sfe)))
  if (length(features_cd)) {
    fd <- colFeatureData(sfe)
    if (!is.null(fd)) {
      out_cd <- .get_not_na_items(fd, features_cd, colname_use)
      features <- setdiff(features, names(out_cd))
    }
  }
  if (!is.null(colGeometryName)) {
    cg <- colGeometry(sfe, colGeometryName, sample_id)
    features_cg <- intersect(features, names(cg))
    if (length(features_cg)) {
      fd <- attr(cg, "featureData")
      if (!is.null(fd)) {
        out_cg <- .get_not_na_items(fd, features_cg, colname_use)
        features <- setdiff(features, names(out_cg))
      }
    }
  }
  if (!is.null(annotGeometryName)) {
    ag <- annotGeometry(sfe, annotGeometryName, sample_id)
    features_ag <- intersect(features, names(ag))
    if (length(features_ag)) {
      fd <- attr(ag, "featureData")
      if (!is.null(fd)) {
        out_ag <- .get_not_na_items(fd, features_ag, colname_use)
        features <- setdiff(features, names(out_ag))
      }
    }
  }
  out <- c(out_rd, out_cd, out_cg, out_ag)
  if (!length(out)) {
    stop("None of the features has the requested metadata.")
  }
  if (length(features)) {
    warning("Features ", paste(features, collapse = ", "),
            " don't have the requested metadata.")
  }
  out
}
