#' @title Patchwork Y-axis
#' 
#' @description Given a table defining a given patchwork contig (generally the 
#' output of patchwork_bins(boundaries_only=TRUE)), define a Y-axis and breaks
#' which are more descriptive of distance between points on the X-axis.
#' 
#' @param tb_region Table defining the regions of the contig, generally the output of patchwork_bins(boundaries_only=TRUE).
#' @param break_num Number of breaks to include on the axis; defaults to 4.
#' 
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import clugPac
#' 
#' @export

patchwork_yaxis   <- function(tb_region,break_num = 4){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange
  
  res <- tb_region %>%
    filter(row_number() == 1) %>%
    mutate(res = (end1 - start1) / bin_alt1_2) %>%
    select(res) %>% 
    unlist(use.names=FALSE)
  brks  <- seq(0,max(tb_region$bin_alt2_2)/2,length.out=break_num) %>% as.integer()
  
  return(setNames(brks,prettyBP(brks * res * 2)))
}