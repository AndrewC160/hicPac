#' @title Rotate pixels 45 degrees.
#'
#' @description Given a table of pixel data, append 'x_coords' and 'y_coords'
#' columns for X and Y values rotated 45 degrees about their center. Useful
#' for plotting HiC data along a horizontal axis without manually rotating
#' plots. Output tables are in long format with duplicated values (each pixel
#' needs four rows, one for each pair of (x,y) coordinates). Adds a "pix_id"
#' column, as well, so ggplot aesthetics should look similar to:
#'
#' ggplot(tb_pix,aes(x=x_coords,y=y_coords,group=pix_id,fill=log10_counts))
#'
#' NOTE: If color/fill values appear to be broken, make sure "group=pix_id" has
#' been specified.
#'
#' While the resulting X-axis should be considerably longer than that of the
#' input plot (2 x resolution^2 = c^2, so pixel width should be
#' sqrt(2 x res^2)), coordinates are divided by sqrt(2) such that rotated pixel
#' diagonals are equal to original pixel sides. This occurs for both X and Y
#' dimensions, and should be consistent for any condition in which matrix is
#' centered on the diagonal.
#'
#' @param df_in Pixel tibble, generally the output format of read_cooler_hdf5().
#'
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import clugPac
#'
#' @export

rotate_pix_45 <- function(df_in) {
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select

  degree  <- pi * -45 / 180
  res     <- (df_in$end_adj1[1] - df_in$start_adj1[1]) / 2
  origin  <- min(df_in$start_adj1)

  #sqrt(2): a^2 + b^2 = c^2;a = b; c = sqrt(2*a^2); for any resolution a ratio is the same: sqrt(2).
  df_in %>%
    #Convert x and y_pos values to doubles to eliminate int64/float64 issues, hopefully...
    mutate_at(c("start_adj1","end_adj1","start_adj2","end_adj2"),as.double) %>%
    mutate(x_pos = (start_adj1 + end_adj1) / 2 - origin,
           y_pos = (start_adj2 + end_adj2) / 2 - origin,
           x_rot = (x_pos * cos(degree) - y_pos * sin(degree)),
           y_rot = (y_pos * cos(degree) + x_pos * sin(degree)),
           x_rot = x_rot / sqrt(2) + origin,
           y_rot = y_rot / sqrt(2) + origin,
           pix_id = paste(bin1_id,bin2_id,sep="_")) %>%
    mutate(x1=x_rot,
           x2=x_rot + res,
           x3=x_rot,
           x4=x_rot - res,
           y1=y_rot + res,
           y2=y_rot,
           y3=y_rot - res,
           y4=y_rot) %>%
    pivot_longer(cols = c("x1","x2","x3","x4","y1","y2","y3","y4"),
                 names_to = c(".value","num"),
                 names_pattern = "([xy])([1234])") %>%
    select(-num,-x_rot,-y_rot) %>%
    rename(x_coords = x,y_coords=y) %>%
    select(-x_pos,-y_pos) %>%
    return()
}
