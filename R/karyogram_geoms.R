#' @title
#' Karyogram geoms
#'
#' @description
#' Return geoms for karyotyping strips on either the X or Y axis.
#'
#' @param gr_in GenomicRanges object representing region to retrive karyotyping data for.
#' @param axis_in Axis to configure geoms for, defaults to "x" and can also be "y".
#' @param range_edge End of the x or y range as appropriate from which to begin drawing geoms.
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import magrittr
#' @import clugPac
#'
#' @export

karyogram_geoms  <- function(gr_in,axis_in="x",range_edge){
  if(axis_in == "y"){
    geom_out  <- get_karyotypes(gr_in) %>%
      as_tibble %>%
      group_split(gieStain) %>%
      lapply(function(tb_y) {
        fill_val  <- tb_y$fill[1]
        geom_rect(data=tb_y,
                  mapping=aes(ymin=adj_start,ymax = adj_end,
                              xmin = range_edge,xmax=Inf),
                  inherit.aes=FALSE,
                  color=NA,fill=fill_val)
      }) %>%
      unlist
  }else{
    geom_out  <- get_karyotypes(grange_in = gr_in) %>%
      as_tibble %>%
      group_split(gieStain) %>%
      lapply(function(tb_x) {
        fill_val  <- tb_x$fill[1]
        geom_rect(data=tb_x,
                  mapping=aes(xmin=adj_start,xmax = adj_end,
                              ymin = range_edge,ymax=Inf),
                  inherit.aes=FALSE,
                  color=NA,fill=fill_val)
      }) %>%
      unlist
  }
  return(geom_out)
}
