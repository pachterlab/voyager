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
    names(out) <- vapply(names(out), .add_sample_id,
        sample_id = sample_id,
        FUN.VALUE = character(1)
    )
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

.add_fd <- function(x, df, res, features, sample_id, name) {
    res <- .add_name_sample_id(res, sample_id)
    df <- .initialize_featureData(df)
    fd <- attr(df, "featureData")
    fd[features, names(res)] <- res
    attr(df, "featureData") <- fd
    df
}

.initialize_fd_dimData <- function(x, MARGIN) {
    fd_name <- switch(MARGIN,
        "rowFeatureData",
        "colFeatureData"
    )
    if (is.null(int_metadata(x)[[fd_name]])) {
        rownames_use <- switch(MARGIN,
            names(rowData(x)),
            names(colData(x))
        )
        fd <- make_zero_col_DFrame(nrow = length(rownames_use))
        rownames(fd) <- rownames_use
        int_metadata(x)[[fd_name]] <- fd
    }
    x
}

# Because adding a new column to S4 DataFrame will remove the attributes
# Put the featureData of colData and rowData in int_metadata instead
.add_fd_dimData <- function(x, MARGIN, res, features, sample_id, type, ...) {
    res <- .add_name_sample_id(res, sample_id)
    x <- .initialize_fd_dimData(x, MARGIN)
    fd_name <- switch(MARGIN,
        "rowFeatureData",
        "colFeatureData"
    )
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

.cbind_all <- function(values, features_list) {
    if (length(values) > 1L) {
        if (length(features_list[["annotgeom"]])) {
            warning("annotGeometry values cannot be cbinded to other values.")
        }
        values <- values[c("assay", "coldata", "colgeom")]
        values <- values[!vapply(values, is.null, FUN.VALUE = logical(1))]
        values <- do.call(cbind, values)
    } else {
        values <- values[[1]]
    }
    values
}

.get_feature_values <- function(sfe, features, sample_id,
                                colGeometryName = NULL, annotGeometryName = NULL,
                                exprs_values = "logcounts", cbind_all = TRUE,
                                show_symbol = TRUE, swap_rownames = "symbol") {
    features_list <- .check_features(sfe, features, colGeometryName, swap_rownames)
    values <- list()
    sample_id_ind <- colData(sfe)$sample_id %in% sample_id
    if (length(features_list[["assay"]])) {
        values_assay <- assay(sfe, exprs_values)[features_list[["assay"]],
            sample_id_ind,
            drop = FALSE
        ]
        # So symbol is shown instead of Ensembl ID
        if (show_symbol && swap_rownames %in% names(rowData(sfe))) {
            rownames(values_assay) <- rowData(sfe)[rownames(values_assay), swap_rownames]
        }
        values_assay <- as.data.frame(as.matrix(t(values_assay)))
        values[["assay"]] <- values_assay
    }
    if (length(features_list[["coldata"]])) {
        values[["coldata"]] <- as.data.frame(colData(sfe)[sample_id_ind,
            features_list[["coldata"]],
            drop = FALSE
        ])
    }
    if (length(features_list[["colgeom"]])) {
        cg <- colGeometry(sfe, colGeometryName, sample_id)
        values[["colgeom"]] <- st_drop_geometry(cg)[sample_id_ind,
            features_list[["colgeom"]],
            drop = FALSE
        ]
    }
    if (length(features_list[["annotgeom"]])) {
        ag <- annotGeometry(sfe, annotGeometryName, sample_id)
        sample_id_ind2 <- ag$sample_id %in% sample_id
        values[["annotgeom"]] <- st_drop_geometry(ag)[sample_id_ind,
            features_list[["annotgeom"]],
            drop = FALSE
        ]
    }
    if (cbind_all) {
        values <- .cbind_all(values, features_list)
    }
    values
}

#' @importFrom SpatialFeatureExperiment localResultFeatures
.check_features_lr <- function(sfe, type, features, sample_id, colGeometryName,
                               annotGeometryName, swap_rownames) {
    features_assay <- localResultFeatures(sfe, type) # includes colData
    features <- .symbol2id(sfe, features, swap_rownames)
    features_assay <- intersect(features_assay, features)
    if (!is.null(colGeometryName)) {
        # Because colGeometryName is also specified for plotting
        features_colgeom <- tryCatch(localResultFeatures(sfe, type,
            colGeometryName = colGeometryName
        ),
        error = function(e) NULL
        )
        features_colgeom <- intersect(features_colgeom, features)
    } else {
        features_colgeom <- NULL
    }
    if (!is.null(annotGeometryName)) {
        features_annotgeom <- tryCatch(localResultFeatures(sfe, type,
            annotGeometryName = annotGeometryName
        ),
        error = function(e) NULL
        )
        features_annotgeom <- intersect(features_annotgeom, features)
    } else {
        features_annotgeom <- NULL
    }
    out <- list(
        assay = features_assay,
        colgeom = features_colgeom,
        annotgeom = features_annotgeom
    )
    other_features <- setdiff(features, Reduce(union, out))
    if (length(other_features)) {
        warning(
            "Features ", paste(other_features, collapse = ", "),
            " are absent in type ", type
        )
    }
    out
}

# Source: https://stackoverflow.com/a/12866609
.unAsIs <- function(X) {
    if("AsIs" %in% class(X)) {
        class(X) <- class(X)[-match("AsIs", class(X))]
    }
    X
}

.get_localResult_attrs <- function(lrs, attribute) {
    out <- lapply(lrs, function(l) {
        if (is.atomic(l) && is.vector(l)) {
            .unAsIs(l)
        } else {
            .unAsIs(l[, attribute])
        }
    })
    data.frame(out, check.names = FALSE)
}

.get_default_attribute <- function(type) {
    switch(type,
        localmoran = "Ii",
        localmoran_perm = "Ii",
        localC_perm = "localC",
        localG = "localG",
        localG_perm = "localG",
        LOSH = "Hi",
        LOSH.mc = "Hi",
        LOSH.cs = "Hi",
        moran.plot = "wx"
    )
}

.get_localResult_values <- function(sfe, type, features, attribute, sample_id,
                                    colGeometryName = NULL,
                                    annotGeometryName = NULL,
                                    cbind_all = TRUE, show_symbol = TRUE,
                                    swap_rownames = "symbol") {
    features_list <- .check_features_lr(
        sfe, type, features, sample_id,
        colGeometryName, annotGeometryName, swap_rownames
    )
    if (is.null(attribute)) {
        attribute <- .get_default_attribute(type)
    }
    values <- list()
    sample_id_ind <- colData(sfe)$sample_id %in% sample_id
    if (length(features_list[["assay"]])) {
        lrs <- localResults(sfe, sample_id, type, features_list[["assay"]])
        if (show_symbol && swap_rownames %in% names(rowData(sfe))) {
            ind <- names(lrs) %in% rownames(sfe)
            names(lrs)[ind] <- rowData(sfe)[names(lrs)[ind], swap_rownames]
        }
        values[["assay"]] <- .get_localResult_attrs(lrs, attribute)
    }
    if (length(features_list[["colgeom"]])) {
        lrs <- localResults(sfe, sample_id, type, features_list[["colgeom"]],
            colGeometryName = colGeometryName
        )
        values[["colgeom"]] <- .get_localResult_attrs(lrs, attribute)
    }
    if (length(features_list[["annotgeom"]])) {
        lrs <- localResults(sfe, sample_id, type, features_list[["annotgeom"]],
            annotGeometryName = annotGeometryName
        )
        values[["annotgeom"]] <- .get_localResult_attrs(lrs, attribute)
    }
    if (cbind_all) {
        values <- .cbind_all(values, features_list)
    }
    values
}

.is_na_list <- function(l) {
    vapply(l, function(d) isTRUE(is.na(d)), logical(1))
}

.get_not_na_items <- function(df, features, colname_use) {
    if (is(df, "sf")) df <- st_drop_geometry(df)
    if (!colname_use %in% names(df)) {
        return(NULL)
    }
    out <- setNames(df[features, colname_use], features)
    out[!.is_na_list(out)]
}

.get_feature_metadata <- function(sfe, features, name, sample_id,
                                  colGeometryName, annotGeometryName,
                                  show_symbol, swap_rownames) {
    colname_use <- .add_sample_id(name, sample_id)
    out_rd <- out_cd <- out_cg <- out_ag <- NULL
    features_rd <- intersect(features, rownames(sfe))
    if (!length(features_rd) && show_symbol && swap_rownames %in% names(rowData(sfe))) {
        features_symbol <- intersect(features, rowData(sfe)[[swap_rownames]])
        .warn_symbol_duplicate(sfe, features_rd)
        features <- setdiff(features, features_symbol)
        features_rd <- rownames(sfe)[match(features_symbol, rowData(sfe)[[swap_rownames]])]
        if (all(is.na(features_rd))) features_rd <- NULL
    }
    if (length(features_rd)) {
        out_rd <- .get_not_na_items(rowData(sfe), features_rd, colname_use)
        features <- setdiff(features, names(out_rd))
        if (show_symbol && swap_rownames %in% names(rowData(sfe))) {
            names(out_rd) <- rowData(sfe)[names(out_rd), swap_rownames]
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
        warning(
            "Features ", paste(features, collapse = ", "),
            " don't have the requested metadata."
        )
    }
    out
}

.local_type2title <- function(type, attribute) {
    if (is.null(attribute)) attribute <- .get_default_attribute(type)
    base <- switch(type,
        localmoran = "Local Moran's I",
        localmoran_perm = "Local Moran's I permutation testing",
        localC = "Local Geary's C",
        localC_perm = "Local Geary's C permutation testing",
        localG = "Getis-Ord Gi(*)",
        localG_perm = "Getis-Ord Gi(*) with permutation testing",
        LOSH = "Local spatial heteroscedasticity",
        LOSH.mc = "Local spatial heteroscedasticity permutation testing",
        LOSH.cs = "Local spatial heteroscedasticity Chi-square test",
        moran.plot = "Moran plot"
    )
    paste0(base, " (", attribute, ")")
}

.deprecate_show_symbol <- function(fun_name, show_symbol, swap_rownames) {
    if (is_present(show_symbol)) {
        deprecate_warn("1.2.0", paste0(fun_name, "(show_symbol = )"),
                       paste0(fun_name, "(swap_rownames = )"))
        # The old behavior
        if (show_symbol) swap_rownames <- "symbol"
    } else show_symbol <- !is.null(swap_rownames)
    list(show_symbol, swap_rownames)
}
