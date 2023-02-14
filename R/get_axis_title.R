#' @title
#' Get axis titles
#'
#' @description
#' Given a tibble returned by read_cooler_hdf5(), axis labels.
#'
#' Typical format:
#' "Chr1: 1,000,000-2,000,000 (1MB, 500kb bin size)" (One chromosome).
#' "300MB, 5MB bin size" (More than one chromosome).
#'
#' @param pix_in Matrix tibble in the format returned by read_matrix_file().
#' @param axis_in Which axis to return labels for. Defaults to "x", can also be "y".
#' @param no_zero_coords Should Zero coordinates be changed to 1? Defaults to TRUE.
#'
#' @import dplyr
#' @import magrittr
#' @import GenomicRanges
#'
#' @export

get_axis_title <- function(pix_in,axis_in = "x",no_zero_coords=TRUE,append_bin_size=TRUE){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select

  if(tolower(axis_in) == "x"){
    tb_pix<- rename(pix_in,
                    seqnames=seqnames1,
                    start = start1,
                    end = end1)
  }else{
    tb_pix<- rename(pix_in,
                    seqnames=seqnames2,
                    start = start2,
                    end = end2)
  }
  tb_pix  <- select(tb_pix,seqnames,start,end) %>%
    distinct %>%
    mutate(width=end - start)
  chroms  <- unique(as.character(tb_pix$seqnames))

  if(append_bin_size){
    bin_size<- paste0("; ",prettyBP(tb_pix$end[1] - tb_pix$start[1])," bin size")
  }else{
    bin_size<- ""
  }

  if(no_zero_coords) tb_pix <- mutate_at(tb_pix, c("start","end"), function(x) ifelse(x == 0,1,x))

  if(length(chroms) == 1){
    gr_out  <- GRanges(chroms,ranges = IRanges(start=min(tb_pix$start),end=max(tb_pix$end)))
    ttl_out <- grange_desc(gr_out,append_ending = bin_size)
  }else{
    wid_val <- tb_pix %>%
      summarize(width = sum(width),.groups="drop") %>%
      unlist(use.names=FALSE) %>%
      prettyBP
    ttl_out <- paste0(wid_val,bin_size)
  }
  return(ttl_out)
}
