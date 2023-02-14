#' @title
#' Plot HiC matrix
#'
#' @description
#' Basic plot of HiC matrix in a region. Accepts a cooler file name and GRanges
#' for X and Y axes, but can also accept a tibble (<tb_pixel>) that has been
#' pre-formatted (generally by a call to read_cooler_hdf5()). In such a case,
#' GRange and <file_cool> arguments are ignored and their values are inferred
#' from the input tibble.
#'
#' @param gr_1 GRange region for X-axis. Defaults to entire genome.
#' @param gr_2 GRange region for Y-axis. Defaults to value of gr_1.
#' @param file_cool Cool file to plot.
#' @param plot_bottom_triangle Should the bottom reflection of the top triangle be included? Defaults to TRUE.
#' @param tb_pixel Tibble formatted from call to read_cooler_hdf5(), can include additional columns but must have basics.
#' @param max_pixel_count Maximum number of pixels to plot, i.e. max_pixels argument for read_cooler_hdf5.
#' @param silent Boolean, defaults to FALSE. Should the function suppress messages about fetching / processing time for read_cooler_hdf5()?
#' @param score_col Column in pixels tibble to plot within the heatmap, defaults to "log10_count".
#' @param fill_low Color of low values (defaults to black).
#' @param fill_high Color of high values (defaults to cyan-ish).
#'
#' @import ggplot2
#' @import tibble
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import clugPac
#' @import GenomicRanges
#'
#' @export

plot_hic_matrix <- function(gr_1=NULL,gr_2=NULL,silent=FALSE,file_cool,tb_pixel=NULL,
                            plot_bottom_triangle=TRUE,max_pixel_count = 6.25e6,
                            score_col = "log10_count",fill_low="#132B43",fill_high="#56B1F7"){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  if(!is.null(tb_pixel)){
    #If pixel table is provided, figure out what the ranges are and go from there.
    gr_1  <- tb_pixel %>%
      group_by(seqnames1) %>%
      summarize(start = min(start1),
                end = max(end1),
                .groups="drop") %>%
      ungroup %>%
      rename(seqnames=seqnames1) %>%
      makeGRangesFromDataFrame(keep.extra.columns = FALSE)
    gr_2  <- tb_pixel %>%
      group_by(seqnames2) %>%
      summarize(start = min(start2),
                end = max(end2),
                .groups="drop") %>%
      ungroup %>%
      rename(seqnames=seqnames2) %>%
      makeGRangesFromDataFrame(keep.extra.columns = FALSE)
    pxls  <- tb_pixel
  }else{
    if(!is.null(gr_1)){
      if(is.null(gr_2)){
        gr_2  <- gr_1
      }
    }
    pxls  <- read_cooler_hdf5(gr_range1 = gr_1,gr_range2=gr_2,file_cool,silent = silent,max_pixels = max_pixel_count)
  }
  pxls    <- mutate(pxls,score = !!as.name(score_col))
  x_rng   <- c(min(pxls$start_adj1),max(pxls$end_adj1))
  y_rng   <- c(min(pxls$start_adj2),max(pxls$end_adj2))
  x_labs  <- get_axis_labs(pix_in = pxls,axis_in = "x")
  y_labs  <- get_axis_labs(pix_in = pxls,axis_in = "y")
  x_ttl   <- get_axis_title(pix_in = pxls,axis_in = "x")
  y_ttl   <- get_axis_title(pix_in = pxls,axis_in = "y",append_bin_size = FALSE)

  p <- ggplot(pxls,aes(xmin=start_adj1,xmax=end_adj1,ymin=start_adj2,ymax=end_adj2,fill=score)) +
    scale_x_continuous(name = x_ttl,expand=c(0,0),limits=x_rng,breaks=x_labs,labels=names(x_labs)) +
    scale_y_continuous(name = y_ttl,expand=c(0,0),limits=y_rng,breaks=y_labs,labels=names(y_labs)) +
    scale_fill_gradient(name = score_col,low = fill_low,high = fill_high,na.value = fill_low) +
    geom_rect(color=NA) +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
          panel.background = element_rect(color=NA,fill=fill_low),
          panel.grid = element_blank(),
          panel.border = element_rect(size=0.5,color="black",fill=NA),
          plot.background = element_rect(color=NA,fill="white"))
  if(plot_bottom_triangle){
    p <- p +
      geom_rect(color=NA,inherit.aes=TRUE,
                mapping=aes(xmin=start_adj2,xmax = end_adj2,
                            ymin=start_adj1,ymax=end_adj1))
  }else{
    p <- p + geom_polygon(
      data=tibble(x=c(x_rng[1],x_rng[2],x_rng[2]),
                  y=c(y_rng[1],y_rng[2],y_rng[1])),
      mapping=aes(x=x,y=y,group=""),inherit.aes = FALSE,
      fill="white",color=NA)
  }
  return(p)
}

