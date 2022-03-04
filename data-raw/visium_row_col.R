# I installed Space Ranger in my home directory
visium_row_col <- read.delim("~/spaceranger-1.3.1/lib/python/cellranger/barcodes/visium-v1_coordinates.txt",
                             header = FALSE, col.names = c("barcode", "col", "row"))
usethis::use_data(visium_row_col, overwrite = TRUE)
