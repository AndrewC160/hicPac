#' @title
#' Plot HiC matrix (pyramid)
#'
#' @description
#' Basic plot of HiC matrix in a region. Accepts a cooler file name and GRanges
#' for X and Y axes, but can also accept a tibble (<tb_pixel>) that has been
#' pre-formatted (generally by a call to read_cooler_hdf5()). In such a case,
#' GRange and <file_cool> arguments are ignored and their values are inferred
#' from the input tibble.
#'
#' NOTE: This function does NOT return standard ggplot objects, it should be
#' aligned/stored using graphics objects like standard grid and gridExtra
#' objects.
#'
#' @param gr_input GRange region for X-axis. Defaults to entire genome.
#' @param file_cool Cool file to plot.
#' @param tb_pixel Tibble formatted from call to read_cooler_hdf5(), can include additional columns but must have basics.
#' @param max_pixel_count Maximum number of pixels to plot, i.e. max_pixels argument for read_cooler_hdf5.
#' @param silent Boolean, defaults to FALSE. Should the function suppress messages about fetching / processing time for read_cooler_hdf5()?
#' @param score_col Column in pixels tibble to plot within the heatmap, defaults to "log10_count".
#' @param score_func Function to apply to score column, for instance < function(x) log10(x+1) >.
#' @param score_name Name to label score with, defaults to the value of <score_col>.
#' @param kary_cm Width of karyogram in fraction of plot height/width; defaults to 2cm.
#' @param as_newpage Should a new page be started for the plot? Defaults to TRUE.
#' @param sub_tracks GGplot with a single panel containing X-axis information.
#' @param sub_track_heigh Ratio of plot area that should be dedicated to the X-axis information panel, defaults to 0.3.
#' @param sub_track_x_ax_frac Fraction of sub_track plot area to dedicate to that plot's x-axis, defaults to 0.2.
#'
#' @import ggplot2
#' @import tibble
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import clugPac
#' @import GenomicRanges
#' @import grid
#' @import gridExtra
#'
#' @export
# cooler_file <- "/N/u/aclugston/projects/hic/hic1/data/384-1/matrix/384-1_5MB.cool"
# tb_pix  <- read_cooler_hdf5(file_cool = cooler_file,silent=FALSE,
#                             cache_tsv = "/N/u/aclugston/projects/eagleC/data/pixels_5MB_384-1.tsv",
#                             max_pixels = Inf,overwrite_cache = FALSE)
# mtx_plot  <- plot_hic_matrix(tb_pixel = tb_pix,score_col = "log10_count",plot_bottom_triangle = FALSE) +
#   scale_shape_manual(values = c("10K"=0,"10K_highres"=1,"50K"=2,"50K_highres"=3,"5K"=4,"5K_combined"=5)) +
#   geom_point(tb_svs,mapping=aes(x=start_adj1,y=start_adj2,shape=res),pch = 0,fill=NA,inherit.aes = FALSE) +
#   theme(panel.border = element_blank())
# kary_plot <- plot_karyotypes(pix_in = tb_pix)
# track_plot<- plot_karyotypes(pix_in = tb_pix,axis_position = "bottom") + theme(axis.title.x = element_text()) + labs(x = "X axis")

plotCombo_pyramid_tracks  <- function(mtx_plot,kary_plot,track_plot=NULL,kary_height=0.05,track_height=0.4,show_viewports=TRUE){
  grid.newpage()
  panel_width <- 0.9
  if(is.null(track_plot)) track_height <- 0
  #After taking into consideration track_height, figure out width of plot (i.e. matrix side length to give
  # a diagonal panel_width wide).
  mtx_side    <- panel_width - track_height
  mtx_diag    <- mtx_side * 2 * 1/sqrt(2)
  mtx_center  <- 1 - mtx_diag/2 - 0.05
  track_height<- mtx_center

  #Having calculated the "actual" matrix center position, adjust track height to reflect it.
  # This should propogate through all dimensions that affect track drawing.
  x_title_frac<- 0.15
  x_text_frac <- 0.15
  x_panel_height<- track_height * (1 - x_title_frac - x_text_frac)
  x_text_height <- track_height * x_text_frac
  x_title_height<- track_height * x_title_frac

  #Define viewports.
  vp_mtx  <- viewport(x=0.5,
                      y=mtx_center,
                      width= unit(mtx_side,"snpc"),
                      height=unit(mtx_side,"snpc"),
                      just=c("center","center"),
                      clip = "off",angle=-45,name = "vp_mtx")
  if(show_viewports) showViewport(vp_mtx)

  #Karyogram.
  vp_kary     <- viewport(
    x=0.5,
    y=unit(mtx_center + mtx_diag/2,"snpc"),
    clip = "off",
    width= unit(mtx_side,"snpc"),
    height=unit(kary_height,"snpc"),
    just=c("right","bottom"),
    angle=45,name="vp_kary")
  if(show_viewports) showViewport(vp_kary)

  #Track panel, background.
  vp_track_back <- viewport(x=0.5,y=mtx_center,
                            width = 1,
                            height = mtx_center,
                            just = c("center","top"),
                            name = "vp_track_bckg")
  if(show_viewports) showViewport(vp_track_back)

  #Track panel.
  vp_track  <- viewport(x=0.5,y=mtx_center,
                      width=unit(mtx_diag,"snpc"),
                      height=unit(x_panel_height,"snpc"),
                      just=c("center","top"),name="vp_track")
  if(show_viewports) showViewport(vp_track)

  #Track x-axis.
  vp_track_x_text <- viewport(x=0.5,
                              y=mtx_center - x_panel_height,
                              width = unit(mtx_diag,"snpc"),
                              height = unit(x_text_height,"npc"),
                              just = c("center","top"),
                              name = "vp_track_x_text")
  if(show_viewports) showViewport(vp_track_x_text)

  #Track x-label.
  vp_track_x_lab  <- viewport(x=0.5,
                              y=mtx_center - x_panel_height - x_text_height,
                              width = unit(mtx_diag,"snpc"),
                              height = unit(x_title_height,"npc"),
                              just = c("center","top"),
                              name = "vp_track_x_lab")
  if(show_viewports) showViewport(vp_track_x_lab)

  #Track y-axis.
  vp_track_y_text  <- viewport(x=unit(0.5,"npc") - unit(mtx_diag,"snpc") * 0.5,
                              y = mtx_center,
                              height = unit(x_panel_height,"npc"),
                              width = unit(0.05,"snpc"),
                              just = c("right","top"),
                              name = "vp_track_y_text")
  if(show_viewports) showViewport(vp_track_y_text)

  #Track y-label.
  vp_track_y_lab  <- viewport(x=(unit(0.5,"npc") - unit(mtx_diag,"snpc") * 0.5) - unit(0.05,"snpc"),
                              y=mtx_center,
                              width = unit(0.05,"snpc"),
                              height = unit(x_panel_height,"npc"),
                              just = c("right","top"),
                              name = "vp_track_y_lab")
  if(show_viewports) showViewport(vp_track_y_lab)

  #Track legend.
  vp_track_legend  <- viewport(x=(unit(0.5,"npc") + unit(mtx_diag,"snpc") * 0.5),
                              y=mtx_center,
                              width = unit(0.1,"snpc"),
                              height = unit(x_panel_height,"npc"),
                              just = c("left","top"),
                              name = "vp_track_legend")
  if(show_viewports) showViewport(vp_track_legend)

  #Legend.
  vp_lgd  <- viewport(x=(unit(0.5,"npc") + unit(mtx_diag,"snpc") * 0.5),
                      y=mtx_center,
                      width=unit(0.1,"npc"),
                      height=unit(0.3,"npc"),
                      just=c("left","bottom"),
                      name="vp_lgd")
  if(show_viewports) showViewport(vp_lgd)

  #Description.
  vp_desc <- viewport(x=0.1,
                      y=0.9,
                      width=mtx_diag/4,
                      height=mtx_diag/4,
                      just=c("left","top"),
                      name = "vp_desc")
  if(show_viewports) showViewport(vp_desc)

  if(!show_viewports){
    #Get and plot sub-track panels/axes
    p_mtx   <- get_panel(mtx_plot)
    p_lgd   <- get_legend(mtx_plot)
    p_desc  <- get_plot_component(p_mtx,"xlab-b")
    pushViewport(vp_mtx)
    grid.draw(p_mtx)
    upViewport()
    pushViewport(vp_lgd)
    grid.draw(p_lgd)
    upViewport()
    pushViewport(vp_desc)
    grid.draw(p_desc)
    upViewport()

    p_kary  <- as_grob(kary_plot)
    pushViewport(vp_kary)
    grid.draw(p_kary)
    upViewport()

    if(!is.null(track_plot)){
      p_track <- get_panel(track_plot)
      p_track_x_txt <- get_plot_component(track_plot,"axis-b")
      p_track_x_ttl <- get_plot_component(track_plot,"xlab-b")
      p_track_y_txt <- get_plot_component(track_plot,"axis-l")
      p_track_y_ttl <- get_plot_component(track_plot,"ylab-l")
      #p_track_legend<- get_plot_component(track_plot,"guide-box")
      p_track_legend<- get_legend(track_plot)

      pushViewport(vp_track_back)
      grid.rect(gp=gpar(col=NA, fill="white"))
      upViewport()
      pushViewport(vp_track)
      grid.draw(p_track)
      upViewport()
      pushViewport(vp_track_x_text)
      grid.draw(p_track_x_txt)
      upViewport()
      pushViewport(vp_track_x_lab)
      grid.draw(p_track_x_ttl)
      upViewport()
      pushViewport(vp_track_y_text)
      grid.draw(vp_track_y_text)
      upViewport()
      pushViewport(vp_track_y_lab)
      grid.draw(p_track_y_ttl)
      upViewport()
      pushViewport(vp_track_legend)
      grid.draw(p_track_legend)
      upViewport()
    }
  }
}

