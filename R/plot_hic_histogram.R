#' @title
#' Plot HiC histogram
#'
#' @description
#' Detect interactions in the genome with a specified ("trap") region and
#' generate a histogram for each bin. Defaults to searching for interactions
#' between the trap region and the rest of its chromosome (<gr_window>="chrom"),
#' but can also search the entire genome (<gr_window>="genome") or a spcific
#' contiguous genomic region (<gr_window>=GRange(...)). Bins are cached to a TSV
#' (<cache_tsv>) rather than the standard pixel file for read_hdf5_*() functions
#' and can be overwritten with <o_write>=TRUE. Function can also bypass the plot
#' output and return a list of tables for generating a similar plot: hist_values
#' contains the histogram, trap describes the trap region, and x_labs contains
#' X-axis labels.
#'
#' @param gr_trap Genomic region to "trap," i.e. what region to detect interactions with.
#' @param gr_window Genomic region within which to search for interactions. Can be a specific GRange, "chrom" (search the trap region's chromosome), or "genome" (search the entire genome).
#' @param cooler_fls Which cooler files to use when determining bins and searching for interactions.
#' @param o_write Should existing <cache_tsv> files be regenerated? Defaults to FALSE.
#' @param retur_plot Should a plot be returned? Defaults to TRUE, and if FALSE a list of tables will be returned for use elsewhere.
#' @param wgs_nrm Should WGS-normalized counts be used? Defaults to TRUE.
#' @param log_scales Should a log10 scale be sued? Defaults to FALSE.
#' @param ignore_diags How many diagonals should be dropped? Defaults to zero (none), while 1 would cause comparisons between the same bin (BinA vs. BinA) to be dropped, 2 would also drop adjacent bins, etc.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import clugPac
#' @import GenomicRanges
#' @import ggplot2
#' @import scales
#' @import GenomicAlignments
#' @import viridis
#'
#' @export

plot_hic_histogram <- function(gr_trap,gr_window="chrom",cooler_fls,o_write=FALSE,return_plot=TRUE,wgs_nrm=TRUE,log_scale=FALSE,max_pix=5E6,ignore_diags=0){
  rename <- dplyr::rename
  mutate <- dplyr::mutate
  filter <- dplyr::filter
  select <- dplyr::select
  seq_szs <- get_seqsizes(as_granges = TRUE)
  if(inherits(gr_window,"character")){
    if(gr_window == "chrom"){
      gr_window <- seq_szs[as.character(seqnames(seq_szs)) == as.character(seqnames(gr_trap))]
    }else
      gr_window <- NULL
  }
  bins<- read_cooler_bins_hdf5(file_cooler = cooler_fls[1]) %>%
    subsetByOverlaps(gr_trap) %>%
    as_tibble %>%
    select(bin_id) %>%
    unlist(use.names=FALSE)

  if(is.null(gr_window)){
    # Genome-wide comparison.
    tb_pix <- read_cooler_hdf5(cooler_fls,max_pixels = max_pix)
  }else{
    # Specific window.
    tb_pix <- read_cooler_hdf5(file_cool = cooler_fls,gr_range1=gr_trap,gr_range2=gr_window,max_pixels = max_pix)
  }
  tb_pix  <- filter(tb_pix,abs(bin1_id-bin2_id) >= ignore_diags)
  if(!"cooler" %in% colnames(tb_pix)) tb_pix <- mutate(tb_pix,cooler="HiC")
  if("count_wgs_nrm" %in% colnames(tb_pix) & wgs_nrm){
    tb_pix <- mutate(tb_pix,count = count_wgs_nrm)
  }
  tb_hist<- tb_pix %>%
    filter(bin1_id %in% bins | bin2_id %in% bins) %>%
    mutate(pre_trap = bin2_id %in% bins) %>%
    mutate(seqnames = ifelse(pre_trap,as.character(seqnames1),as.character(seqnames2)),
           start = ifelse(pre_trap,start1,start2),
           end = ifelse(pre_trap,end1,end2),
           start_adj = ifelse(pre_trap,start_adj1,start_adj2),
           end_adj = ifelse(pre_trap,end_adj1,end_adj2)) %>%
    select(seqnames,start,end,start_adj,end_adj,count,bin1_id,bin2_id,cooler) %>%
    group_by(seqnames,start,end,start_adj,end_adj,cooler) %>%
    reframe(count = sum(count),
            bin1_id = ifelse(bin1_id > bin2_id,bin1_id,bin2_id),
            bin2_id = ifelse(bin1_id > bin2_id,bin2_id,bin1_id))

  seq_adj <- get_seqsizes_adj()[as.character(seqnames(gr_trap))]
  tb_trap <- tb_hist %>%
    group_by(cooler) %>%
    summarize(max_count=1.1*max(count),.groups="drop") %>%
    mutate(start = start(gr_trap),
           end = end(gr_trap),
           start_adj=start + seq_adj,
           end_adj = end + seq_adj,
           mid = (start + end) / 2,
           mid_adj = (start_adj + end_adj)/2) %>%
    select(cooler,start,end,mid,start_adj,end_adj,mid_adj,max_count)

  bin_size<- prettyBP(tb_hist$end[1] - tb_hist$start[1])
  if(length(unique(tb_hist$seqnames)) > 1){
    tb_labs <- tb_hist %>%
      group_by(seqnames) %>%
      summarize(start = min(start),
                end=max(end),
                start_adj = min(start_adj),
                end_adj = max(end_adj),
                .groups="drop") %>%
      arrange(start_adj) %>%
      mutate(shade = row_number() %% 2 == 0)

    x_labs  <- vectify(tb_labs,start_adj,seqnames)
    scale_x <- list(scale_x_continuous(name=paste0("Genomic bins (",bin_size," bins)"),
                                       expand=c(0,0),breaks=x_labs,labels=names(x_labs),
                                       oob=scales::oob_keep),
                    geom_rect(data=filter(tb_labs,shade),ymin=-Inf,ymax=Inf,fill="blue",alpha=0.1))
  }else{
    x_labs  <- tb_hist %>%
      mutate(ax_bin = cut(start_adj,breaks=4,labels = FALSE)) %>%
      ungroup %>%
      group_by(ax_bin) %>%
      filter(row_number() == 1) %>%
      ungroup %>%
      mutate(label=prettyBP(start)) %>%
      vectify(start_adj,label)
    scale_x <- scale_x_continuous(name=grange_desc(gr_window,append_ending = paste(";",bin_size,"bin size")),
                                  expand=c(0,0),limits=c(min(tb_hist$start_adj),max(tb_hist$end_adj)),
                                  breaks=x_labs,labels=names(x_labs),oob=scales::oob_keep)
  }

  if(!return_plot){
    out_p <- list(hist_values=tb_hist,
                  trap = tb_trap,
                  labels = x_labs)
  }else{
    if(log_scale){
      scale_y <- scale_y_log10(name="Interactions",expand=c(0,0),labels=function(x) prettyNumbers(x,drop_zero = TRUE))
    }else{
      scale_y <- scale_y_continuous(name="Interactions",expand=c(0,0),labels=function(x) prettyNumbers(x,drop_zero = TRUE))
    }

    out_p <- ggplot(tb_hist,aes(xmin=start_adj,xmax=end_adj,ymin=0,ymax=count,fill=log10(count+1))) +
      facet_wrap(.~cooler,ncol=1,scales="free_y",strip.position = "right") +
      scale_x +
      scale_y +
      scale_fill_viridis(option="magma") +
      geom_rect(color=NA) +
      geom_rect(data=tb_trap,mapping=aes(xmin=start_adj,xmax=end_adj,ymin=-Inf,ymax=Inf),color="red",fill=NA,linewidth=0.5,inherit.aes=FALSE) +
      geom_point(data=tb_trap,mapping=aes(x=mid_adj,y=1.1*max_count),color="red",pch=18,size=5,inherit.aes=FALSE) +
      theme(plot.background = element_rect(fill="white",color=NA),
            panel.background = element_blank(),
            panel.border = element_rect(fill=NA,color="black",linewidth=0.5),
            panel.grid = element_blank(),
            panel.spacing.y = unit(0,"lines"),
            strip.background = element_blank(),
            strip.text.y.right = element_text(angle=0),
            axis.text.x = element_text(angle=45,hjust=1,vjust=1),
            axis.ticks.x = element_blank(),
            legend.position = "none")
  }
  return(out_p)
}
