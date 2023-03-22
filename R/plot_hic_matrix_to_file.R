#' @title
#' Plot HiC matrix to file.
#'
#' @description
#' Executes plot_hic_matrix(), but does so with a file as output. Useful for
#' caching plot outputs in markdowns, and has the added benefit of configuring
#' asymmetric matrices them with the correct aspect ratios. To specify
#' approximate figure size, either height or width can be specified. If both
#' are provided, these will override the aspect ratio calculation and use these
#' dimensions directly. Output returns a file name. Note: height/width values
#' apply to the panel itself, not plot margins or axes. As such, a figure with
#' <fig_width> of 6 inches will be saved as ~7.5-8 inches wide depending on
#' labels and legends.
#'
#' @param gr_1 GRange region for X-axis. Defaults to entire genome.
#' @param gr_2 GRange region for Y-axis. Defaults to value of gr_1.
#' @param file_cool Cool file to plot.
#' @param tb_pixel Tibble formatted from call to read_cooler_hdf5(), can include additional columns but must have basics.
#' @param kary_width Width of karyogram in fraction of plot height/width; defaults to 0.25 (2.5\% of the total plot height and width are added as a karyogram).
#' @param max_pixel_count Maximum number of pixels to plot, i.e. max_pixels argument for read_cooler_hdf5.
#' @param silent Boolean, defaults to FALSE. Should the function suppress messages about fetching / processing time for read_cooler_hdf5()?
#' @param score_col Column in pixels tibble to plot within the heatmap, defaults to "log10_count".
#' @param file_output Output file name. If NULL, no file will be cached.
#' @param fig_width Desired width of file. If not provided, will default to 6 (inches implied).
#' @param fig_height Desired height of file.
#' @param fig_units Unit input for ggsave, defaults to "inches".
#' @param fig_device Device input for ggsave, defaults to "png".
#' @param fig_dpi DPI input for ggsave, defauts to 300.
#' @param overwrite_cache Should the function ignore existing cached images and overwrite them? Defaults to FALSE.
#'
#' @import ggplot2
#' @import egg
#' @import grid
#' @import tibble
#' @import dplyr
#'
#' @export

plot_hic_matrix_to_file <- function(gr_1=NULL,gr_2=NULL,silent=FALSE,file_cool=NULL,tb_pixel=NULL,kary_width=0.025,max_pixel_count = 6.25e6,score_col = "log10_count",
                            file_output=NULL,fig_width=NULL,fig_height=NULL,fig_units="in",fig_device="png",fig_dpi=300,overwrite_cache=FALSE){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter

  if(is.null(file_output)){
    stop("No output file provided, exiting.")
  }

  if(!file.exists(file_output) | overwrite_cache){
    fig_dims  <- get_figure_dims(gr_1=gr_1,gr_2=gr_2,fig_height=fig_height,fig_width = fig_width)
    plot_hic_matrix(gr_1=gr_1,gr_2=gr_2,silent=silent,file_cool=file_cool,tb_pixel=tb_pixel,kary_width=kary_width,max_pixel_count = max_pixel_count,score_col=score_col) %>%
      set_panel_size(width=unit(fig_dims$width,fig_units),height = unit(fig_dims$height,fig_units)) %>%
      ggsave(filename=file_output,device = fig_device,dpi = fig_dpi,height= 1.3*fig_dims$height,width=1.3*fig_dims$width)
  }
  return(file_output)
}

