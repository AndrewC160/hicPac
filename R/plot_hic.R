#' @title
#' Plot HiC
#'
#' @description
#' Plot a standard HiC interaction matrix (input format should be that returned
#' by read_hic_matrix()).
#'
#' @param matrix_tbl Tibble returned by read_hic_matrix().
#' @param reflect_xy Draw both (i,j) and (j,i) tiles (reflect matrix across diagonal i=j). Defaults to FALSE.
#'
#' @import dplyr
#' @import ggplot2
#' @import clugPac
#'
#' @export

plot_hic  <- function(matrix_tbl,reflect_xy=FALSE){
  mtx_tb  <- matrix_tbl
  x_labs  <- get_axis_labs(mtx_in = mtx_tb,axis_in = "x") %>%
    vectify(plot_bin,label)
  y_labs  <- get_axis_labs(mtx_tb,axis_in = "y") %>%
    vectify(plot_bin,label)
  #For whole-genome tables, instead of representing (i,j) and (j,i) both in a table 2x longer,
  #just include both an (i,j) and a (j,i) tile geom.
  if(reflect_xy){
    geom_reflect  <- geom_tile(aes(x=binY,y=binX))
  }else{
    geom_reflect  <- NULL
  }
  x_rng   <- range(mtx_tb$binX)
  y_rng   <- range(mtx_tb$binY)

  ggplot(mtx_tb,aes(x=binX,y=binY,fill=log10(count))) +
    scale_x_continuous(expand = expansion(mult=c(0,0)),
                       breaks = x_labs,
                       labels = names(x_labs)) +
    scale_y_continuous(expand = expansion(mult=c(0,0)),
                       breaks = y_labs,
                       labels = names(y_labs)) +
    scale_fill_gradient() +
    geom_tile() +
    geom_reflect +
    theme(panel.background = element_rect(fill="black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank())
}
