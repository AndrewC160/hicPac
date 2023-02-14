#' @title
#' Read Cooler HDF5, rotate 45 degrees.
#'
#' @description
#' Read .cool files as HDF5 using the read_cooler_hdf5() function, but perform
#' required steps to rotate its pixels 45 degrees within a given window. In
#' short: 1) Read cooler file(s) with read_cooler_hdf5(), 2) rotate pixels with
#' rotate_pix_45(), 3) trim pixels to the specified X range with
#' trim_rot45().
#'
#' @param file_cool Name of .cool file(s) with the appropriate bin size. If multiple files are provided, functions is run recursively and tables are combined (including cache TSV).
#' @param gr_range1 GRange of contiguous region for X dimension. Defaults to entire genome, and must be contiguous.
#' @param diag_distace Maximum distance from the diagonal (useful for filtering pixels from tables meant for rotated plots). Filter is applied AFTER fetching pixels, so does not generally help performance.
#' @param silent Should processing time messages be suppressed? Defaults to FALSE.
#' @param max_pixels Maximum number of pixels to retrieve without erroring out. Defaults to 6.25 million (2500x2500 grid), can also be set to Inf to disregard.
#' @param cache_table Filename of TSV to cache pixel results into. If provided and this file exists, this table will be read and returned. If not found, it will be created.
#' @param overwrite_cache Should cache file be overwritten? Defaults to FALSE.
#'
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import data.table
#' @import GenomicRanges
#' @import magrittr
#' @import tictoc
#' @import clugPac
#' @import rhdf5
#' @import Rhdf5lib
#' @import data.table
#'
#' @export

read_cooler_hdf5_rot45  <- function(file_cool,gr_range1=NULL,diag_distance=NULL,silent=TRUE,max_pixels=6.25e6,cache_tsv=NULL,overwrite_cache = FALSE){
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  rename  <- dplyr::rename

  read_cooler_hdf5(file_cool=file_cool,
                   gr_range1 = gr_range1,
                   diag_distance = diag_distance,
                   silent=silent,
                   max_pixels = max_pixels,
                   cache_tsv=cache_tsv,
                   overwrite_cache=overwrite_cache) %>%
    rotate_pix_45 %>%
    trim_rot45(gr_x = gr_range1) %>%
    return()
}
