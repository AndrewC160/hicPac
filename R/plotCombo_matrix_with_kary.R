#' @title
#' Arrangement of square matrix panel with karyograms on top and right axes.
#'
#' @description
#' Provide a matrix plot and two karyogram plots (one for the X and Y axes),
#' returns a combined plot with the karyograms inserted.
#'
#' @param p_mtx GGplot object for matrix. Axis labels and titles will be extracted from this plot.
#' @param p_xkary GGplot object for X axis karyogram.
#' @param p_ykary GGplot object for Y axis karyogram.
#' @import ggplot2
#' @import clugPac
#' @import cowplot
#' @import grid
#' @import gridExtra
#'
#' @export

plotCombo_matrix_with_kary <- function(p_mtx,p_xkary,p_ykary){
  mtx_panel <- get_panel(p_mtx)
  mtx_xlab  <- get_plot_component(p_mtx,"xlab-b")
  mtx_xaxis <- get_plot_component(p_mtx,"axis-b")
  mtx_ylab  <- get_plot_component(p_mtx,"ylab-l")
  mtx_yaxis <- get_plot_component(p_mtx,"axis-l")
  mtx_legend<- get_plot_component(p_mtx,"guide-box")
  pkx       <- get_panel(p_xkary)
  pky       <- get_panel(p_ykary)

  grid.arrange(
    pkx,       #1
    mtx_ylab,  #2
    mtx_yaxis, #3
    mtx_panel, #4
    pky,       #5
    mtx_legend,#6
    mtx_xaxis, #7
    mtx_xlab,  #8
    layout_matrix = as.matrix(
      rbind(
        c(NA,NA,1,NA,NA,NA),
        c(2,3,4,5,NA,6),
        c(NA,NA,7,NA,NA,NA),
        c(NA,NA,8,NA,NA,NA)
      )
    ),
    heights=c(1,30,2,1),
    widths =c(1,1,30,1,1,3)
  )
}

