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
#' @importFrom SingleCellExperiment int_colData int_colData<-
#' @importFrom methods new
#' @importFrom utils packageVersion
.add_name_sample_id <- function(out, sample_id) {
    names(out) <- vapply(names(out), .add_sample_id,
        sample_id = sample_id,
        FUN.VALUE = character(1)
    )
    out
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
        if (is.character(type)) type2 <- get(type, mode = "S4")
        attribute <- info(type2, "default_attr")
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
                                  reducedDimName,
                                  show_symbol, swap_rownames) {
    colname_use <- .add_sample_id(name, sample_id)
    out_rd <- out_cd <- out_cg <- out_ag <- out_reddim <- NULL
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
    if (!is.null(reducedDimName)) {
        reddim <- reducedDim(sfe, reducedDimName)
        if (is.null(colnames(reddim))) {
            names_use <- paste0(reducedDimName, ncol(reddim))
        } else names_use <- colnames(reddim)
        if (is.numeric(features))
            features_reddim <- names_use[features]
        else
            features_reddim <- intersect(features, names_use)
        fd <- reducedDimFeatureData(sfe, reducedDimName)
        if (!is.null(fd)) {
            out_reddim <- .get_not_na_items(fd, features_reddim, colname_use)
            if (is.numeric(features)) features <- NULL
            else
                features <- setdiff(features, names(out_reddim))
        }
    }
    out <- c(out_rd, out_cd, out_cg, out_ag, out_reddim)
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
    if (is.character(type)) type <- get(type, mode = "S4")
    if (is.null(attribute)) attribute <- info(type, "default_attr")
    base <- info(type, "title")
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

#' @importFrom SingleCellExperiment reducedDim reducedDim<-
.add_reddim_colnames <- function(sfe, dimred) {
    rd <- reducedDim(sfe, dimred)
    if (is.null(colnames(rd)))
        colnames(rd) <- paste0(dimred, seq_len(ncol(rd)))
    reducedDim(sfe, dimred) <- rd
    sfe
}

# As in MatrixExtra, only for Csparse for now
.empty_dgc <- function(nrow, ncol) {
    out <- new("dgCMatrix")
    out@Dim <- as.integer(c(nrow, ncol))
    out@p <- integer(ncol+1L)
    out
}
