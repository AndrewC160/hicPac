#' @title
#' Plot karyotypes
#'
#' @description
#' Plot karyotyping facet for a given genomic range. If no range is provided,
#' draw for entire genome.
#'
#' @param grange_win GRanges object to illustrate. Assumes contiguous regions; otherwise use mutliple plot calls.
#' @param pix_in Pixel tibble in the format of read_cooler_hdf5() output. If provided, this will be used to determine <grange_win>.
#' @param pix_axis Axis of pixel tibble for which to plot karyogram, defaults to "x" (i.e. seqnames1,start1,start2), can also be "y".
#' @param on_axis Axis on which to plot karyogram, defaults to "x" and can also be "y".
#' @param break_num Number of breaks to split axis into. If multiple chromosomes are included, the start of each chromosome will also be included, so this value is approximate.
#' @param axis_position Position of axis; defaults to "top" in which case x-axis is drawn on the top of the panel. If axis is drawn in y dimension, this is changed to "left".
#'
#' @import ggplot2
#' @import GenomicRanges
#' @import dplyr
#' @import tibble
#' @import magrittr
#'
#' @export

plot_karyotypes <- function(grange_win=NULL,pix_in=NULL,pix_axis="x",on_axis="x",break_num = 5,axis_position="top"){
  seq_adjs  <- get_seqsizes_adj()
  if(!is.null(pix_in)){
    if(tolower(pix_axis) == "x"){
      tb_pix  <- rename(pix_in,seqnames=seqnames1,start=start1,end=end1)
    }else{
      tb_pix  <- rename(pix_in,seqnames=seqnames2,start=start2,end=end2)
    }
    grange_win <- tb_pix %>%
      select(seqnames,start,end) %>%
      distinct %>%
      mutate(start_adj = start + seq_adjs[as.character(seqnames)],
             end_adj = end + seq_adjs[as.character(seqnames)]) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  }

  if(is.null(grange_win)){
    grange_win <- get_seqsizes(as_granges=TRUE)
  }
  grange_win$start_adj  <- start(grange_win) + seq_adjs[as.character(seqnames(grange_win))]
  grange_win$end_adj    <- end(grange_win) + seq_adjs[as.character(seqnames(grange_win))]

  x_rng   <- c(min(grange_win$start_adj),max(grange_win$end_adj))

  tb_bands  <- get_karyotypes() %>%
    as_tibble %>%
    filter(start_adj >= x_rng[1],
           end_adj <= x_rng[2])

  seq_num  <- length(unique(tb_bands$seqnames))
  if(seq_num > break_num){
    ax_labs<- tb_bands %>%
      group_by(seqnames) %>%
      summarize(start_adj = min(start_adj),
                end_adj = max(end_adj),
                start = min(start),
                .groups="drop") %>%
      mutate(label = ifelse(start == 1,
                            as.character(seqnames),
                            paste0(as.character(seqnames),":",
                              prettyNum(start,big.mark = ",")))) %>%
      vectify(start_adj,label)
  }else{
    ax_bins <- seq(min(tb_bands$start_adj),max(tb_bands$end_adj),length.out=break_num + 1)
    ax_labs <- tb_bands %>%
      mutate(bin = cut(start_adj,breaks = ax_bins,include.lowest=TRUE)) %>%
      group_by(seqnames,bin) %>%
      summarize(start = min(start),
                end = max(end),
                start_adj = min(start_adj),
                end_adj = max(end_adj),
                .groups="drop") %>%
      group_by(seqnames) %>%
      mutate(first_seq = row_number() == 1) %>%
      ungroup %>%
      mutate(label = case_when(start == 1 ~ as.character(seqnames),
                               first_seq ~ paste0(as.character(seqnames),prettyNum(start,big.mark=",")),
                               TRUE ~ prettyNum(start,big.mark=","))) %>%
      vectify(start_adj,label)
  }

  if(tolower(on_axis) == "x"){
    p <- ggplot(tb_bands,aes(xmin=start_adj,xmax=end_adj,ymin=0,ymax=1,fill=fill)) +
      scale_x_continuous(limits=x_rng,expand=c(0,0),breaks = ax_labs,labels = names(ax_labs),position = axis_position) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_identity() +
      geom_rect(color=NA) +
      theme(legend.position = "none",
            axis.text.x.top = element_text(angle=-45,hjust=0.5,vjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(size=0.5,color="black",fill=NA),
            plot.margin = unit(c(0,0,0,0),units = "lines"))
  }else{
    if(axis_position == "top") axis_position  <- "left"
    p <- ggplot(tb_bands,aes(ymin=start_adj,ymax=end_adj,xmin=0,xmax=1,fill=fill)) +
      scale_y_continuous(limits=x_rng,expand=c(0,0),breaks = ax_labs,labels = names(ax_labs),position = axis_position) +
      scale_x_continuous(expand=c(0,0)) +
      scale_fill_identity() +
      geom_rect(color=NA) +
      theme(legend.position = "none",
            axis.text.y = element_text(hjust=1,vjust=0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            strip.text.x = element_text(angle=0,hjust=0,vjust=0.5),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(size=0.5,color="black",fill=NA),
            plot.margin = unit(c(0,0,0,0),units = "lines"))
  }
  return(p)
}
