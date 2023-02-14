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

plot_hic_matrix_pyramid <- function(gr_input=NULL,silent=FALSE,file_cool,tb_pixel=NULL,as_newpage=TRUE,
                                    sub_tracks = NULL,sub_track_height=0.3,sub_track_x_ax_frac = 0.1,
                                    max_pixel_count = 6.25e6,score_col = "log10_count",score_func=NULL,score_name=NULL,kary_cm=2){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  if(!is.null(tb_pixel)){
    #If pixel table is provided, figure out what the ranges are and go from there.
    gr_input  <- tb_pixel %>%
      group_by(seqnames1) %>%
      summarize(start = min(start1),
                end = max(end1),
                .groups="drop") %>%
      ungroup %>%
      rename(seqnames=seqnames1) %>%
      makeGRangesFromDataFrame(keep.extra.columns = FALSE)
    pxls  <- tb_pixel
  }else{
    pxls  <- read_cooler_hdf5(gr_range1 = gr_input,file_cool,silent = silent,max_pixels = max_pixel_count)
  }
  pxls  <- mutate(pxls,score = !!as.name(score_col))

  if(!is.null(score_func)){
    pxls<- mutate(pxls,score = score_func(score))
  }

  score_nm<- ifelse(is.null(score_name),score_col,score_name)

  x_labs  <- get_axis_labs(pix_in = pxls,axis_in = "x")

  bin_size<- paste0("; ",prettyBP(pxls$end_adj1[1] - pxls$start_adj1[1])," bin size")
  x_ttl   <- ifelse(is.null(gr_input),"",grange_desc(gr_input,append_ending = bin_size))

  x_rng   <- c(min(pxls$start_adj1),max(pxls$end_adj1))

  p_kary  <- get_karyotypes(grange_in = gr_input) %>%
    as_tibble %>%
    ggplot(aes(xmin=adj_start,xmax=adj_end,ymin=0,ymax=1,fill=fill)) +
      scale_x_continuous(name = x_ttl,limits=x_rng,expand=c(0,0),position="top",labels=names(x_labs),breaks=x_labs) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_identity() +
      geom_rect(color=NA,size=0.25) +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(size=0.25,color="black",fill=NA),
            axis.text.x = element_text(angle=-45,hjust=1,vjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = unit(c(0,0,0,-0.001),"npc"))

  p_mtx   <- ggplot(pxls,aes(xmin=start_adj1,xmax=end_adj1,ymin=start_adj2,ymax=end_adj2,fill=score)) +
    scale_x_continuous(name = x_ttl,expand=c(0,0),limits=x_rng,breaks=x_labs,labels=names(x_labs)) +
    scale_y_continuous(expand=c(0,0),limits=x_rng) +
    scale_fill_gradient(name = score_col,low="white",high="red") +
    geom_rect(color=NA) +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(color=NA,fill="white"))
  g_mtx   <- get_panel(p_mtx)
  g_lgd   <- get_legend(p_mtx)

  #Grid viewports (matrix, legend, and karyogram).
  if(is.null(sub_tracks)) sub_track_height  <- 0
  x_ax_height   <- sub_track_x_ax_frac * sub_track_height
  x_pan_height  <- sub_track_height - x_ax_height

  mtx_side_len  <- (0.85 - sub_track_height) * 2 * 1/sqrt(2)
  trk_side_len  <- sqrt(2*mtx_side_len^2)
  mtx_side_len  <- unit(mtx_side_len,"snpc")
  trk_side_len  <- unit(trk_side_len,"snpc")

  vp_mtx  <- viewport(x=0.5,y=sub_track_height,
                      width= mtx_side_len,
                      height=mtx_side_len,
                      just=c("center","center"),
                      clip = "off",angle=-45,name = "matrix.vp")
  vp_kar_h<- ifelse(x_ttl == "",kary_cm,1.3*kary_cm) %>% unit("cm")
  vp_karL <- viewport(x=0.5,y=0.85,
                      clip = "off",
                      width= mtx_side_len,
                      height=vp_kar_h,
                      just=c("right","bottom"),
                      angle=45,name="karyogram_left.vp")
  vp_lgd  <- viewport(x=1,y=0.9,width=unit(0.2,"snpc"),height=unit(0.3,"snpc"),just=c("right","top"),name="legend.vp")
  vp_plt  <- viewport(x=0.5,y=sub_track_height,
                      width=trk_side_len,
                      height=unit(x_pan_height,"snpc"),
                      just=c("center","top"),name="tracks.vp")
  vp_plt_y_left   <- viewport(x=unit(0.5,"npc") - trk_side_len * 0.5,
                              y=sub_track_height,
                              width=unit(3,"cm"),
                              height=unit(x_pan_height,"snpc"),
                              just = c("right","top"),
                              name = "tracks_yaxis_left.vp")
  vp_plt_y_right  <- viewport(x=unit(0.5,"npc") + trk_side_len * 0.5,
                              y=sub_track_height,
                              width=unit(3,"cm"),
                              height=unit(x_pan_height,"snpc"),
                              just = c("left","top"),
                              name = "tracks_yaxis_right.vp")
  vp_plt_x_bottom <- viewport(x=unit(0.5,"npc"),
                              y=0,
                              width = trk_side_len,
                              height=unit(x_ax_height,"snpc"),
                              just = c("center","bottom"),
                              name = "tracks_xaxis_bottom.vp")
  grid.newpage()
  showViewport(vp_mtx)
  showViewport(vp_lgd)
  showViewport(vp_karL)
  showViewport(vp_plt)
  showViewport(vp_plt_y_left)
  showViewport(vp_plt_y_right)
  showViewport(vp_plt_x_bottom)

  #Get sub-track panels/axes
  g_strack_panel    <- get_panel(sub_tracks)
  g_strack_ax_left  <- cowplot::get_y_axis(sub_tracks,position = "left")
  g_strack_ax_right <- cowplot::get_y_axis(sub_tracks,position = "right")
  g_strack_ax_bot   <- cowplot::get_x_axis(sub_tracks,position = "bottom")
  g_strack_lgd      <- cowplot::get_legend(sub_tracks)

  # Draw plot components in their respective viewports.
  suppressWarnings({
    if(as_newpage) grid.newpage()
    pushViewport(vp_mtx)
    grid.draw(g_mtx)
    upViewport()
    pushViewport(vp_lgd)
    grid.draw(g_lgd)
    upViewport()
    print(p_kary,vp=vp_karL)
    pushViewport(vp_plt)
    grid.draw(g_strack_panel)
    upViewport()
    pushViewport(vp_plt_y_left)
    grid.draw(g_strack_ax_left)
    upViewport()
    pushViewport(vp_plt_y_right)
    grid.draw(g_strack_ax_right)
    upViewport()
    pushViewport(vp_plt_x_bottom)
    grid.draw(g_strack_ax_bot)
    upViewport()
  })
}

