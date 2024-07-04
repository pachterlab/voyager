.get_tx_df <- function(sfe, data_dir, tech, file, sample_id, spatialCoordsNames,
                       gene_col, bbox, gene, return_sf = FALSE,
                       rowGeometryName = "txSpots", geoparquet_file = NULL,
                       flip = FALSE) {
    if (is.null(sfe) && is.null(data_dir) && is.null(file)) {
        stop("One of sfe, data_dir, and file must be specified.")
    }
    if (is.null(sfe)) {
        if (!is.null(data_dir)) {
            c(spatialCoordsNames, gene_col, cell_col, file) %<-%
                getTechTxFields(tech, data_dir)
        }
        else if (!is.null(file)) file <- normalizePath(file, mustWork = TRUE)
        if (grepl("\\.parquet$", file)) {
            check_installed("arrow")
            df <- arrow::read_parquet(file)
        } else {
            check_installed("data.table")
            df <- data.table::fread(file)
        }
        if (!identical(gene, "all")) {
            gene <- gene[gene %in% df[[gene_col]]]
            if (!length(gene)) stop("None of the genes are in the transcript spot file. Please check gene_col.")
            df <- df[df[[gene_col]] %in% gene,]
        }
        df <- df[,c(spatialCoordsNames[1:2], gene_col)]
        names(df) <- c("X", "Y", "gene")
        if (flip) df$Y <- -df$Y
    } else {
        if (!identical(gene, "all")) {
            gene <- gene[gene %in% rownames(sfe)]
            if (!length(gene)) stop("None of the genes are in the SFE object.")
        }
        spatialCoordNames <- c("X", "Y")
        sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
        if (!is.null(geoparquet_file)) {
            if (!gdalParquetAvailable())
                stop("GDAL Parquet driver is required to selectively read genes from GeoParquet file.")
            if (length(sample_id) > 1L) {
                if (!any(sample_id %in% names(geoparquet_file)))
                    stop("GeoParquet file name vector must have names to match them to the sample_ids.")
            } else {
                if (length(geoparquet_file) > 1L) {
                    geoparquet_file <- geoparquet_file[1]
                    message("Only the first GeoParquet file specified will be used.")
                }
                names(geoparquet_file) <- sample_id
            }
        }
        dfs <- lapply(sample_id, function(s) {
            rgn <- .check_rg(rowGeometryName, sfe, s)
            sql_req <- gdalParquetAvailable() && !is.na(geoparquet_file[s]) &&
                !identical(gene, "all")
            if (rgn %in% rowGeometryNames(sfe)) {
                tx <- rowGeometry(sfe, type = rowGeometryName, sample_id = s)
                if (!identical(gene, "all")) tx <- tx[rownames(sfe) %in% gene,]
                # always "gene" if read from formatTxSpots
                which_empty <- st_is_empty(tx)
                if (any(which_empty) && sql_req) {
                    do_sql <- TRUE
                    gene_sql <- rownames(tx)[which_empty]
                    tx <- tx[!which_empty,]
                } else do_sql <- FALSE
            } else if (sql_req) {
                do_sql <- TRUE
                gene_sql <- gene
                tx <- NULL
            } else {
                stop(rowGeometryName, " not found in the SFE object")
            }
            if (do_sql) {
                geoparquet_file <- normalizePath(geoparquet_file, mustWork = TRUE)
                tx_sql <- readSelectTx(geoparquet_file, gene_sql)
                rownames(tx_sql) <- tx_sql$gene
                tx <- rbind(tx, tx_sql)
            }
            if (return_sf) {
                out <- tx
            } else {
                out <- as.data.frame(st_coordinates(tx))
            }
            out$sample_id <- s
            out
        })
        df <- do.call(rbind, dfs) # data.table indexing geometry takes a while
    }
    df <- .crop(df, bbox)
    df
}

#' Plot transcript spot density as 2D histogram
#'
#' @inheritParams plotCellBin2D
#' @inheritParams plotGeometry
#' @inheritParams SpatialFeatureExperiment::getTechTxFields
#' @inheritParams SpatialFeatureExperiment::aggregateTx
#' @param data_dir Top level directory of the output files. This can be
#' specified in place of \code{sfe} to directly read the transcript spot
#' coordinates from the file. When reading from file, transcripts from this file
#' are plotted and argument \code{sample_id} is ignored.
#' @param file File (not GeoParquet) with numeric columns for xy coordinates of
#' the transcript spots. Ignored if \code{data_dir} is specified.
#' @param gene Character vector of names of genes to plot. If "all" then
#' transcript spots of all genes are plotted.
#' @param flip Logical, whether to flip the y axis when plotting data from file.
#' @return A ggplot object, facetting by sample.
#' @export
#' @importFrom zeallot %<-%
#' @importFrom SpatialFeatureExperiment .check_rg
plotTxBin2D <- function(sfe = NULL, data_dir = NULL,
                        tech = c("Vizgen", "Xenium", "CosMX"), file = NULL,
                        sample_id = "all", bins = 200, binwidth = NULL,
                        hex = FALSE, ncol = NULL, bbox = NULL, gene = "all",
                        spatialCoordsNames = c("X", "Y"), gene_col = "gene",
                        rowGeometryName = "txSpots", flip = FALSE, tx_file = NULL) {
    if (!is.null(sfe)) {
        sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
        tech <- match.arg(tech)
    }
    if (length(sample_id) > 1L && !is.null(sfe)) use_samples <- TRUE else use_samples <- FALSE
    df <- .get_tx_df(sfe, data_dir, tech, file, sample_id, spatialCoordsNames,
                     gene_col, bbox, gene, rowGeometryName = rowGeometryName,
                     flip = flip, geoparquet_file = tx_file)
    p <- ggplot(df, aes(X, Y))
    if (hex) p <- p + geom_hex(bins = bins, binwidth = binwidth)
    else p <- p + geom_bin2d(bins = bins, binwidth = binwidth)
    p <- p + scale_fill_distiller(palette = "Blues", direction = 1) +
        coord_equal() +
        scale_x_continuous(expand = expansion()) +
        scale_y_continuous(expand = expansion()) +
        labs(x = NULL, y = NULL)
    if (use_samples) p <- p + facet_wrap(~ sample_id, ncol = ncol)
    p
}
