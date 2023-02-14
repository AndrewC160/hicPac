#' @title Rotate points 45 degrees.
#' 
#' @description Given a matrix with x and y columns, and a range describing an
#' X and y axis, rotate points in the input matrix 45 degrees and rescale them
#' such that the new x-axis (originally the diagonal) is the same length as the
#' original x-axis (originally the side).
#' 
#' @param tb_in Tibble input.
#' @param x_grange GRanges object (one contiguous range) describing the limits of the x-axis.
#' @param y_grange GRanges object (one contiguous range) describing the limits of the y-axis.
#' @param x_rng Vector of length 2 with max and min values of x-range.
#' @param y_rng Vector of length 2 with max and min values of y-range.
#' @param x_col Column name of x-values.
#' @param y_col Column name of y-values.
#' @param x_col_new Name of new rotate column. Defaults to NULL, in which case x_col is replaced with the new value.
#' @param y_col_new Name of new rotate column. Defaults to NULL, in which case y_col is replaced with the new value.
#' 
#' import tibble
#' 
#' @export

rotate_points_45<- function(tb_in,x_grange,y_grange,x_rng=NULL,y_rng=NULL,x_col=NULL,y_col=NULL,x_col_new = NULL,y_col_new = NULL){
  rename <- dplyr::rename
  if(is.null(x_col)) stop("No x input column name provided.")
  if(is.null(y_col)) stop("No y input column name provided.")
  if(!is.null(x_grange)){
    x_rng <- c(start(x_grange),end(x_grange))
  }
  if(!is.null(y_grange)){
    y_rng <- c(start(y_grange),end(y_grange))
  }
  if(is.null(x_rng)) stop("No x-range provided.")
  if(is.null(y_rng)) stop("No y-range provided.")
  
  degree  <- pi * -45 / 180
  origin  <- c(x_rng[1],y_rng[1])
  x_diff  <- diff(x_rng)
  y_diff  <- diff(y_rng)
  
  #Steps: 1) Define points on scale of 0-1 on each axis.
  #       2) Rotate about artificial origin of (0,0).
  #       3) Divide by sqrt(2) to resize diagonal to axis size.
  #       4) Multiply fraction of new rotated axis and add origin back.
  tb_out <- tb_in %>%
    mutate(x_frac = (!!as.name(x_col) - origin[1]) / x_diff,
           y_frac = (!!as.name(y_col) - origin[2]) / y_diff,
           x_rot = (x_frac * cos(degree) - y_frac * sin(degree))/sqrt(2),
           y_rot = (y_frac * cos(degree) + x_frac * sin(degree))/sqrt(2),
           x_rot = x_rot * x_diff + origin[1],
           y_rot = y_rot * y_diff + origin[2]) %>%
    select(-x_frac,-y_frac)
  
  if(is.null(x_col_new) | is.null(y_col_new)){
    tb_out <- mutate(tb_out,
                     !!as.name(x_col) := x_rot,
                     !!as.name(y_col) := y_rot) %>%
      select(-x_rot,-y_rot)
  }else{
    tb_out <- rename(tb_out,
                     !!as.name(x_col_new) := x_rot,
                     !!as.name(y_col_new) := y_rot)
  }
  return(tb_out)
}
