#' @title
#' Patchwork rotate 45 degrees
#' 
#' @description 
#' Given a matrix of patchwork pixels (typically the output of 
#' patchwork_matrix()), convert all bins to pixels and rotate their coordinates
#' by 45 degrees about the left-most pixels in the contig in question. Produces
#' a table 4x longer than its input (one line per corner of each pixel), and 
#' produces <x_coords>, <y_coords>, and <pix_id> columns that can be used to 
#' plot in GGPlot2, i.e. ggplot(mtx,aes(x=x_coords,y=y_coords,group=pix_id)).
#' 
#' @param mtx_in Patchwork matrix produced by patchwork_matrix() to be rotated.
#' 
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' 
#' @export

patchwork_rot45 <- function(mtx_in){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange
  
  degree  <- pi * -45 / 180
  
  #sqrt(2): a^2 + b^2 = c^2;a = b; c = sqrt(2*a^2); for any resolution a ratio is the same: sqrt(2).
  mtx_in %>%
    mutate(x_pos = bin_alt1 - 0.5, #Start at the center of the pixel.
           y_pos = bin_alt2 - 0.5,
           x_rot = (x_pos * cos(degree) - y_pos * sin(degree)),
           y_rot = (y_pos * cos(degree) + x_pos * sin(degree)),
           x_rot = x_rot / sqrt(2),#Scale the new coordinates back to the original dimensions (i.e. shrink rotated triangle).
           y_rot = y_rot / sqrt(2),
           pix_id = paste(bin_alt1,bin_alt2,sep="_")) %>%
    mutate(x1=x_rot, #Identify each coordinate for pixel polygons.
           x2=x_rot + 0.5*sqrt(2), #0.5*sqrt(2) = half of diagonal for unit pixels.
           x3=x_rot,
           x4=x_rot - 0.5*sqrt(2),
           y1=y_rot + 0.5*sqrt(2),
           y2=y_rot,
           y3=y_rot - 0.5*sqrt(2),
           y4=y_rot) %>%
    pivot_longer(cols = c("x1","x2","x3","x4","y1","y2","y3","y4"),
                 names_to = c(".value","num"),
                 names_pattern = "([xy])([1234])") %>%
    mutate_at(c("x","y"),round,digits=3) %>%
    mutate_at(c("x","y"),function(x) ifelse(x < 0,0,x)) %>%
    select(-num,-x_rot,-y_rot,-x_pos,-y_pos) %>%
    rename(x_coords = x,y_coords=y) %>%
    return()
}