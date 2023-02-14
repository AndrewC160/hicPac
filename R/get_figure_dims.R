#' @title 
#' Get figuire dimensions
#' 
#' @description 
#' Given a pair of genomic ranges indicating X and Y axes fore a matrix, return
#' a vector with appropriate height and width values that reflect the true
#' dimensions. If one or both GRanges are NULL, the aspect ratio is set to 1.
#' In the event that one axis is more than <max_difference> times the length of
#' the other, a warning is displayed and dimensions are returned for the 
#' <max_difference> ratio. Can be overridden by increasing <max_difference>.
#' 
#' @param gr_1 GRange region for X-axis. Defaults to entire genome.
#' @param gr_2 GRange region for Y-axis. Defaults to value of gr_1.
#' @param fig_height Desired height of figure (units don't matter).
#' @param fig_width Desired width of figure (units don't matter).
#' @param max_difference Maximum ratio difference between height and width.
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export

get_figure_dims  <- function(gr_1=NULL,gr_2=NULL,fig_height=NULL,fig_width=NULL,max_difference=20){
  if(is.null(gr_1) | is.null(gr_2)){
    asp_ratio <- 1
  }else{
    asp_ratio   <- width(gr_1) / width(gr_2)
    if(asp_ratio > max_difference){
      message("Aspect ratio of ",round(asp_ratio,digits=2)," greater than max_difference, defaulting to ",max_difference,".")
      asp_ratio <- max_difference
    }else if(asp_ratio < 1/max_difference){
      message("Aspect ratio of ",round(asp_ratio,digits=2)," less than  1/max_difference, defaulting to ",round(1/max_difference,digits=2),".")
    }
  }
  if(is.null(fig_height) & is.null(fig_width)){
    fig_width   <- 6
    fig_height  <- fig_width / asp_ratio
  }else if(is.null(fig_height)){
    fig_height  <- fig_width / asp_ratio
  }else if(is.null(fig_width)){
    fig_width   <- fig_height * asp_ratio
  }
  
  return(list(height=fig_height,width=fig_width))
}