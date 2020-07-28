#' S4 class for SpatialImage with transcript locations and cell segmentations
#'
#' The SMSeurat class inherits from the \code{SpatialImage} class from Seurat, and
#' extends it to smFISH based or ISS based datasets where transcript locations
#' and/or cell segmentation polygons are available.
#'
#' @importClassesFrom Seurat SpatialImage
#' @slot coordinates A data frame with 3 columns: "x", and "y" for
#' the x and y coordinates, "cell" for cell IDs.
#' @slot tissue_owin A list of one \code{spatstat::owin} object, or the
#' observation window of the tissue.
#' @slot qhulls A list whose length is the same as the number of rows of
#' coordinates, for cell segmentation polygons.
#' @slot spots A list whose length is the same as the number of rows of
#' coordinates, for transcript coordinates. Each element is a ppp.
#' @slot unitname Name of the unit of the coordinates.
#' @export
#' @importFrom spatstat owin
SMSeurat <- setClass(
  Class = "SMSeurat",
  contains = "SpatialImage",
  slots = list(
    coordinates = "data.frame",
    tissue_owin = "list",
    qhulls = "list",
    spots = "list",
    unitname = "character"
  )
)

setValidity("SMSeurat", function(object) {
  if (!is.data.frame(object@coordinates)) {
    return("Slot coordinates must be a data frame.")
  } else if (!is.character(object@unitname) || length(object@unitname) > 1) {
    return("Slot unitname must be a character vector of length 1.")
  } else if (!is.list(object@qhulls)) {
    return("Slot qhulls must be a list.")
  } else if (!all(c("x", "y", "cell") %in% names(object@coordinates))) {
    return("Slot coordinates must have columns x, y, and cell")
  } else if (!is.list(object@spots)) {
    return("Slot spots must be a list.")
  } else if (length(object@tissue_owin) > 1 || !is(object@tissue_owin[[1]], "owin")) {
    return("Slot tissue_owin's only element must be an owin object.")
  } else if (length(object@qhulls) != nrow(object@coordinates)) {
    return("Length of qhulls must be the same as the number of rows of coordinates.")
  } else if (length(object@spots) != nrow(object@coordinates)) {
    return("Length of spots must be the same as the number of rows of coordinates.")
  } else {
    cls <- map(object@qhulls, class)
    cls <- Reduce(intersect, cls)
    if (!"data.frame" %in% cls) {
      return("Each entry of qhulls must be a data frame.")
    }
  }

  nms <- map(object@qhulls, names)
  nms <- Reduce(intersect, nms)
  if (!all(c("x", "y") %in% nms)) {
    return("For each cell, the length of the x coordinate vector of qhulls must be the
           same as the length of the y coordinate vector.")
  }
  cls <- unique(map_chr(object@spots, class))
  if (cls != "ppp") {
    return("Each element of spots must be a ppp object.")
  }
  TRUE
})

#' Create a SMSeurat object
#'
#' This function takes in data frames and converts the transcript locations
#' into ppp. The SMSeurat object is used for data from smFISH or ISS based methods
#' that has transcript coordinates and better yet cell segmentation polygons.
#' This can be used to analyze subcellular transcript localization and use the
#' subcellular localization to generate cell level features.
#'
#' @param cell_centroids A data frame with 3 columns. The first column is
#' assumed to be the x coordinate, and the second column is assumed to be the y
#' coordinate. The third column is cell IDs. This is not turned into a ppp
#' object for ease of plotting with ggplot2, since I foresee that I'll plot cell
#' locations more often than transcript spots. However, this will likely change
#' over the development of this package after I write functions to plot various
#' operations on ppp objects that only have base R plotting methods in spatstat.
#' @param tissue_owin An \code{spatstat::owin} object, or the observation window
#' of the tissue, such as the extent of a field of view or the boundary of the
#' tissue section on the slide. This is optional as this is not always provided
#' by published data. If omitted, and if cell_segmentations provided, then this
#' will be the convex hull of all cell segmentations. If cell_segmentations is
#' also omitted, then this will be the convex hull of all the transcript coordinates.
#' @param tx_coords A data frame with at least these 4 columnns: The first column
#' is assumed to be the x coordinate of the transcript spot, the second one is
#' assumed to be the y coordinate, the third one the cell this transcript
#' belongs to, and the fourth the gene the transcript is assigned to. All other
#' columns will become marks (i.e. metadata or covariate) of the points. Again,
#' @param cell_segmentations A data frame with 3 columns: The first is assumed
#' to be the x coordinate of the segmentation boundnaries, the second is the y
#' coordinate, and the third the cell ID. This is optional since for published
#' datasets, this is not always available. When cell segmentations are unavailable,
#' the convex hull of the transcript points are used. See the spatstat book for
#' why this is problematic. The points should go counterclockwise and the first
#' vertex should not be repeated.
#' @param assay Character, the assay of the Seurat object this spatial info
#' corresponds to.
#' @param key Character, must end with "_", the key Seurat use for plots.
#' @param unitname Unit to use
#' @importFrom dplyr group_nest mutate rename semi_join
#' @importFrom magrittr %>%
#' @importFrom purrr map map2
#' @importFrom grDevices chull
#' @importFrom methods new
#' @importFrom spatstat ppp
#' @importFrom tidyr unnest
#' @export
CreateSMSeurat <- function(cell_centroids, tx_coords, tissue_owin = NULL,
                           cell_segmentations = NULL, unitname = "um",
                           assay = "Spatial", key = "spatial_") {
  names(cell_centroids)[1:3] <- c("x", "y", "cell")
  names(tx_coords)[1:4] <- c("x", "y", "cell", "gene")
  mark_nms <- setdiff(names(tx_coords), c("x", "y", "cell"))
  unitname <- as.character(unitname[1])
  tx_coords <- tx_coords %>%
    group_nest(cell, .key = "coords") %>%
    semi_join(cell_centroids, by = "cell")
  if (nrow(tx_coords) < 1L) {
    stop("Cells in tx_coords do not match those in cell_centroids.")
  }
  if (nrow(tx_coords) < nrow(cell_centroids)) {
    n <- nrow(cell_centroids) - nrow(tx_coords)
    cell_centroids <- cell_centroids %>%
      semi_join(tx_coords, by = "cell")
    warning(n, " cells do not have transcript locations; these cells are removed.\n")
  }
  tx_coords <- tx_coords[match(cell_centroids$cell, tx_coords$cell),]
  if (!is.null(cell_segmentations)) {
    names(cell_segmentations)[1:3] <- c("x", "y", "cell")
    cell_segmentations <- cell_segmentations %>%
      group_nest(cell, .key = "coords")
  } else {
    # Find convex hull, problematic, don't use unless you really don't
    # have any info about it from published data.
    cell_segmentations <- tx_coords %>%
      # chull gives coordinates in clockwise order but owin requires counterclockwise
      mutate(hull_coords = map(coords, function(df) df[rev(chull(df)),])) %>%
      select(-coords) %>%
      rename(coords = hull_coords)
  }
  cell_segmentations <- cell_segmentations %>%
    semi_join(cell_centroids, by = "cell")
  if (nrow(cell_segmentations) < 1L) {
    stop("Cells in cell_segmentations do not match those in cell_centroids.")
  }
  if (nrow(cell_segmentations) < nrow(cell_centroids)) {
    n <- nrow(cell_centroids) - nrow(cell_segmentations)
    cell_centroids <- cell_centroids %>%
      semi_join(cell_segmentations, by = "cell")
    tx_coords <- tx_coords %>%
      semi_join(cell_segmentations, by = "cell")
    warning(n, "cells do not have cell segmentations; these cells and their transcripts are removed.\n")
  }
  cell_segmentations <- cell_segmentations[match(cell_centroids$cell,
                                                 cell_segmentations$cell),]
  if (is.null(tissue_owin)) {
    tissue_df <- cell_segmentations %>%
      unnest(cols = "coords")
    inds <- rev(chull(tissue_df[, c("x", "y")]))
    tissue_coords <- tissue_df[inds, c("x", "y")]
    tissue_owin <- owin(poly = tissue_coords)
  }
  # Construct the pattern column
  spots <- map2(tx_coords$coords, cell_segmentations$coords,
                function(.x, .y) {
                  ppp(x = .x$x, y = .x$y, poly = .y, marks = .x[, mark_nms])
                })
  new("SMSeurat", coordinates = cell_centroids,
      qhulls = as.list(cell_segmentations$coords),
      tissue_owin = list(tissue_owin), spots = as.list(spots),
      unitname = unitname, assay = assay, key = key)
}

#' Set coordinates of an SM object
#'
#' @param object An SM object.
#' @importFrom tibble tibble
#' @export
SetCoords <- function(object, new_coords) {
  object@coordinates <- new_coords
  object
}
# To do: Getters and setters of all slots of SM
#' Method for SM
#'
#' @param object The SM object.
#' @param rownames Logical, whether cell IDs should be rownames or a column on
#' its own. I hate rownames of data frames; no wonder setting rownames of tibbles
#' has been deprecated.
#' @export
GetTissueCoordinates.SMSeurat <- function(object, qhulls = FALSE, rownames = TRUE) {
  if (!qhulls) {
    out <- object@coordinates
    if (rownames) {
      out <- as.data.frame(out)
      rownames(out) <- out$cell
      out$cell <- NULL
    }
  } else {
    out <- tibble(cell = object@coordinates$cell,
                  coords = object@qhulls) %>%
      unnest(cols = "coords")
  }
  out
}

#' Get the transcript locations
#'
#' The title is self-explanatory
#'
#' @param object An SMSeurat object.
#' @export
setGeneric("GetSpots", function(object) standardGeneric("GetSpots"))

#' @rdname GetSpots
#' @export
setMethod("GetSpots", signature = "SMSeurat",
          function(object) object@spots)

#' Add image to Seurat object
#'
#' I haven't found a function in Seurat that does this.
#'
#' @param object An object.
#' @param image An object that inherits from the class SpatialImage.
#' @export
#' @importFrom methods is validObject
#' @importFrom Seurat Cells

setGeneric("AddImage",
           function(object, image) standardgGeneric(object, image),
           signature = c("object", "image"))

#' @rdname AddImage
setMethod("AddImage", signature = c("Seurat", "SMSeurat"),
          function(object, image) {
            if (any(!image@assay %in% Seurat::Assays(object))) {
              stop("Image must be associated with an existing assay in the Seurat object.")
            }
            # It's Seurat team's business to implement getter functions
            img_name <- gsub("_$", "", image@key[1])
            # Check that the coordinates match the number and IDs of cells
            cds <- GetTissueCoordinates(image, rownames = FALSE)
            if (nrow(cds) != ncol(object)) {
              stop("Number of cell coordinate pairs must be the same as the number of cells.")
            } else if (any(!cds$cell %in% Seurat::Cells(object))) {
              stop("Cell IDs in image must match those in the Seurat object.")
            } else {
              cds <- cds[match(Seurat::Cells(object), cds$cell),]
              image <- SetCoords(image, cds)
              validObject(image)
            }
            object@images[[img_name]] <- image
            object
          })

#' @method dim SMSeurat
#' @export
dim.SMSeurat <- function(x) {
  coords <- GetTissueCoordinates(object = x)
  return(c(
    max(coords[, 1]) - min(coords[, 1]),
    max(coords[, 2]) - min(coords[, 2])
  ))
}

#' @method GetImage SMSeurat
#' @export
GetImage.SMSeurat <- function(
  object,
  mode = c('grob', 'raster', 'plotly', 'raw'),
  ...
) {
  mode <- match.arg(arg = mode)
  return(NullImage(mode = mode))
}

#' @method Cells SMSeurat
#' @export
Cells.SMSeurat <- function(x) {
  coords <- GetTissueCoordinates(x, rownames = FALSE)
  coords$cell
}

#' @method RenameCells SMSeurat
#' @export
RenameCells.SMSeurat <- function(x, new.names = NULL) {
  if (!is.null(new.names)) {
    stopifnot(length(new.names) == nrow(x@coordinates))
  }
  x@coordinates$cell <- new.names
  x
}

#' @method Radius SMSeurat
#' @export
#'
Radius.SMSeurat <- function(object) {
  return(NULL)
}

#' @method subset SMSeurat
#' @export
#'
subset.SMSeurat <- function(x, cells, ...) {
  coordinates <- GetTissueCoordinates(object = x, rownames = FALSE)
  coordinates <- coordinates %>%
    filter(cell %in% cells)
  slot(object = x, name = 'coordinates') <- coordinates
  return(x)
}

#' @export
#' @method [ SMSeurat
"[.SMSeurat" <- function(x, i, ...) {
  return(subset(x = x, cells = i))
}

# Return a null image
#
# Copies straight from Seurat. I don't like this; this means that the spatial
# part of Seurat is not very extensible. And this is the release CRAN version!
# @param mode Image representation to return
# see \code{\link{GetImage}} for more details
#
#' @importFrom grid nullGrob
#' @importFrom grDevices as.raster
#
NullImage <- function(mode) {
  image <- switch(
    EXPR = mode,
    'grob' = nullGrob(),
    'raster' = as.raster(x = new(Class = 'matrix')),
    'plotly' = list('visible' = FALSE),
    'raw' = NULL,
    stop("Unknown image mode: ", mode, call. = FALSE)
  )
  return(image)
}

