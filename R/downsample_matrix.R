#' @title
#' Downsample HiC matrix.
#'
#' @description
#' Given a HiC matrix (as read by read_matrix_file()), downsample it to have
#' approximately <bins_X> by <bins_Y> bins. If there fewer than or equal to
#' that many bins in the original matrix, return the original matrix.
#'
#' @import dplyr
#' @import tibble
#' @import data.table
#'
#' @export

downsample_matrix <- function(matrix_in,x_bins = 750,y_bins = 750,return_extra_cols=FALSE){
  #matrix_in = hic_main()
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange

  #Figure out per-chromosome bin allotment.
  x_brks <- matrix_in %>%
    group_by(seqnamesX) %>%
    summarize(width = max(endX) - min(startX),.groups="drop") %>%
    mutate(frac = width / sum(width),
           num_bins = ceiling(x_bins * frac),
           lag_bins = cumsum(lag(num_bins,default = 0)))
  x_lags <- vectify(x_brks,lag_bins,seqnamesX)
  x_brks <- vectify(x_brks,num_bins,seqnamesX)

  y_brks <- matrix_in %>%
    group_by(seqnamesY) %>%
    summarize(width = max(endY) - min(startY),.groups="drop") %>%
    mutate(frac = width / sum(width),
           num_bins = ceiling(y_bins * frac),
           lag_bins = cumsum(lag(num_bins,default = 0)))
  y_lags <- vectify(y_brks,lag_bins,seqnamesY)
  y_brks <- vectify(y_brks,num_bins,seqnamesY)

  #Sort into X and Y bins, respectively, and summarize.
  mtx_out <- matrix_in %>%
    group_by(seqnamesX) %>%
    mutate(plot_binX = cut(binX,breaks = x_brks[as.character(first(seqnamesX))],labels=FALSE,include.lowest=TRUE)) %>%
    ungroup %>%
    group_by(seqnamesY) %>%
    mutate(plot_binY = cut(binY,breaks = y_brks[as.character(first(seqnamesY))],labels=FALSE,include.lowest=TRUE)) %>%
    ungroup %>%
    group_by(seqnamesX,seqnamesY,plot_binX,plot_binY) %>%
    summarize(startX = min(startX),
              startY = min(startY),
              endX = max(endX),
              endY = max(endY),
              start_binX = min(binX),
              end_binX = max(binX),
              start_binY = min(binY),
              end_binY = max(binY),
              count = mean(count),
              .groups="drop") %>%
    mutate(binX = plot_binX + x_lags[as.character(seqnamesX)],
           binY = plot_binY + y_lags[as.character(seqnamesY)]) %>%
    select(-plot_binX,-plot_binY)
  if(return_extra_cols){
    mtx_out   <- select(mtx_out,seqnamesX,startX,endX,seqnamesY,startY,endY,binX,binY,count,everything())
  }else{
    mtx_out   <- select(mtx_out,seqnamesX,startX,endX,seqnamesY,startY,endY,binX,binY,count)
  }
  return(mtx_out)
}
