#' @title Distance-normalize pixels.
#'
#' @description
#' Given a (non-rotated) pixel table from read_cooler_hdf5(), add columns for
#' distance-normalized count values (counts to normalize can be provided via
#' <count_col>, but should be non-logarithmic values, see citation.). Function
#' will operate in groups based on chromosomes and samples (coolers), and will
#' separate distances into logarithmic bins (see The Hitchhiker’s guide to Hi-C
#' analysis: Practical guidelines, Methods, Brian Lajoie, et.al., 2015) which
#' double in size in series (1MB, 2MB, 4MB, etc.).
#'
#' Function can return only the break values (bins_only=TRUE), or can return a
#' tibble with distance-normalized values. Default behavior returns the input
#' pixel table (tb_pix_in) with new columns "dist" (distance in BP between bins,
#' NA if bins are on different chromosome), "dist_bin" (the bin in which the
#' pixel falls), "dist_count" (count column with distance baseline subtracted
#' and negative values zapped to zero), "dist_base" (baseline of interactions
#' among pixels in the same distance bin).
#'
#' @param tb_pix_in Square matrix pixel table from read_cooler_hdf5() (not rotated).
#' @param count_col Column of count values to normalize. Defaults to "count".
#' @param bins_only Boolean; should only a vector of bin breaks be returned? Defaults to FALSE.
#' @param table_only Boolean; should only a table summarizing breaks and baseline values per cell line and chromosome be returned? Defaults to FALSE; superseded by bins_only.
#'
#' @import dplyr
#' @import magrittr
#' @import clugPac
#' @import tidyr
#'
#' @export

distance_normalize_pixels <- function(tb_pix_in,count_col = "count",bins_only=FALSE,table_only=FALSE){
  bin_num   <- 8
  dist_bins <- 1E6 * 2^c(0:bin_num) - 1E6
  dist_bins <- setNames(dist_bins[-(bin_num+1)],
                        paste0(prettyBP(dist_bins[-(bin_num+1)]),"-",
                               prettyBP(dist_bins[-1])))
  if(bins_only){
    dat_out <- dist_bins
  }else{
    if(!"cooler" %in% colnames(tb_pix_in)){
      tb_pix_in <- mutate(tb_pix_in,"cooler")
    }
    tb_p  <- tb_pix_in %>%
      mutate(dist = ifelse(seqnames1 != seqnames2,NA,(end2 + start2)/2 - (end1 + start1)/2)) %>%
      mutate(dist_bin = cut(dist,right = TRUE,breaks = dist_bins,include.lowest = TRUE,labels = FALSE),
             dist_bin = names(dist_bins)[dist_bin],
             dist_bin = ifelse(is.na(dist_bin),"Trans",dist_bin),
             dist_bin = factor(dist_bin,levels=c(names(dist_bins),"Trans")),
             dist_count = !!as.name(count_col))

    tb_base <- tb_p %>%
      group_by(cooler,seqnames1,seqnames2,dist_bin) %>%
      summarize(dist_base = mean(dist_count),.groups="drop")

    if(table_only){
      dat_out <- tb_base
    }else{
      dat_out <- left_join(tb_p,tb_base,by=c("cooler","seqnames1","seqnames2","dist_bin")) %>%
        mutate(dist_count = dist_count - dist_base) %>%
        mutate(dist_count = ifelse(dist_count < 0,0,dist_count))
    }
  }
  return(dat_out)
}


