#' My alternative to Seurat's plotting function
#'
#' Because at present, Seurat's version is not extensible, though at present,
#' my version is simpler. For now it only colors by the metadata. It would also
#' only work on SMSeurat. Or maybe I can contribute some spatial Seurat plotting
#' alternatives to dittoSeq.
#'
#' @param object A Seurat object
#' @param group.by Same as in Seurat's SpatialPlot.
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot geom_polygon aes
#' @importFrom rlang !! sym
#' @export
SpatialPlot2 <- function(object, group.by = NULL) {
  image_use <- object[[Images(object)[1]]]
  cell_segs <- GetTissueCoordinates(image_use, qhulls = TRUE)
  # There should be a hierarchy in the SpatialImage class inheritance
  # One level down is smFISH or ISS vs. array.
  # Because of my AddImage, they should be in the same order.
  metas <- object@meta.data %>%
    mutate(cell = rownames(.))
  cell_segs <- cell_segs %>%
    left_join(metas, by = "cell")
  if (length(group.by) > 1L) {
    # To do; now just do 1 feature
    stop("More than one variables is not supported yet.")
  } else {
    ggplot(cell_segs, aes(x, y, fill = !!sym(group.by), group = cell)) +
      geom_polygon() +
      coord_equal()
  }
}

SpatialFeaturePlot2 <- function(object, features) {
  image_use <- object[[Images(object)[1]]]
  df <- GetAssayData(object)[features,] %>%
    as.matrix() %>%
    as.data.frame()
  df$cell <- Cells(object)
  cell_segs <- GetTissueCoordinates(image_use, qhulls = TRUE)
  cell_segs <- cell_segs %>%
    left_join(df, by = "cell")
  if (length(features) > 1L) {
    # To do; now just do 1 feature
    stop("More than one variables is not supported yet.")
  } else {
    ggplot(cell_segs, aes(x, y, fill = V1, group = cell)) +
      geom_polygon() +
      coord_equal()
  }
}
