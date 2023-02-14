#' @title Trim rotated pixels
#'
#' @description Given a tibble with rotated pixels (output of rotate_pix_45()),
#' perform two functions. First, if a grange object has been provided as
#' <gr_x>, remove all rows where the x-axis values are outside of this range.
#' Second, remove all rows which would fall outside of the x and y axes in
#' in general. This has the effect of "shaving" top,bottom,and side "corners"
#' from each pixel such that the top/bottom/side boundaries are flat, i.e. the
#' pixels become triangles. Looks cleaner in plots, plus GGplot's clipping
#' system doesn't like to draw partial geoms. If no arguments are provided
#' besides the input tibble, edge pixels are determined and trimmed
#' appropriately.
#'
#' @param tb_in Tibble input, generally output of rotate_pix_45().
#' @param gr_x If provided, all pixels will be filtered to keep those which fall inside this genomic region.
#' @param x_range Vector of upper and lower bounds to use manually on x axis, i.e. the "start_adj1/end_adj1" boundaries to filter by. Ignored if gr_x is provided.
#' @param y_range Vector of upper and lower bounds to use manually on y axis, i.e. the "start_adj2/end_adj2" boundaries to filter by.
#'
#' @import tibble
#' @import tidyr
#' @import data.table
#' @import magrittr
#' @import GenomicRanges
#' @import clugPac
#'
#' @export

trim_rot45  <- function(tb_in,gr_x=NULL,x_range=NULL,y_range=NULL){
  between   <- data.table::between
  filter    <- dplyr::filter
  if(!is.null(gr_x)){
    sq_nm   <- as.character(seqnames(gr_x))
    seq_adj <- get_seqsizes_adj()[sq_nm]
    x_range <- c(seq_adj + start(gr_x),seq_adj + end(gr_x))
  }
  if(is.null(x_range)) x_range  <- range(tb_in$x_coords)
  if(is.null(y_range)) y_range  <- range(tb_in$y_coords)

  px_diag <- sqrt(2 * (tb_in$end1[1] - tb_in$start1[1])^2)

  tb_in %>%
    filter(between(y_coords,y_range[1],y_range[2]-px_diag/4,incbounds = FALSE),
           between(x_coords,x_range[1]+px_diag/4,x_range[2]-px_diag/4,incbounds = FALSE)) %>%
    return()
}
