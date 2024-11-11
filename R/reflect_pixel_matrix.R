#' @title Reflect pixel matrix.
#' 
#' @description
#' Reflect a (non-rotated) pixel matrix across the diagonal (i.e. swap columns
#' with start/end/seqnames values). Note that the only columns that are swapped
#' are the seqnames, start, end, start_adj, end_adj, and bin1/bin2 id columns.
#' To use rotate_pix_45(), perform the reflection first, i.e.:
#' 
#' read_cooler_hdf5(<file_in) %>%
#'  reflect_pixel_matrix() %>%
#'  rotate_pix_45()
#' 
#' @param tb_pix_in Tibble of pixel data (as output by read_cooler_hdf5(), typically).
#' 
#' @import tibble
#' @import clugPac
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' 
#' @export

reflect_pixel_matrix  <- function(tb_pix_in){
  tb_pix_in %>%
    swap_columns("seqnames1","seqnames2") %>%
    swap_columns("start1","start2") %>%
    swap_columns("end1","end2") %>%
    swap_columns("start_adj1","start_adj2") %>%
    swap_columns("end_adj1","end_adj2") %>%
    swap_columns("bin1_id","bin2_id") %>%
    return()
}
