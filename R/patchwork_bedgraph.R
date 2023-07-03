#' @title Patchwork BedGraph
#'
#' @description
#' Plot signal from BedGraph files within a genome patchwork, respecting
#' segment orientation. Returns a GGplot object.
#'
#' @param regions_tb Region table, generally the output of patchwork_bins(boundaries_only=TRUE).
#' @param bdg_files Filenames of Tabix-indexed BedGraphs. Assumes index file is in place.
#' @param bdg_weights Weights to be multiplied against scores in each bedGraph file. Defaults to 1.
#' @param bdg_names Names for each BedGraph file, generally their signal contents (i.e. epitope).
#' @param epitope_colors Specify colors for each included epitope to be passed to ggplot2's scale_fill_manual(); i.e. (K27Ac="green",K4Me3="blue"). Should be named similar to bdg_files, but doesn't need to be.
#' @param bdg_units Y-scale title; defaults to "Pileup."
#' @param raw_alpha Alpha (transparency) between 0 and 1 for raw scores; defaults to 0.5.
#' @param weighted_alpha Alpha (transparency) between 0 and 1 for weight-multiplied scores; defaults to 1.
#' @param num_y_brks Number of Y-axis breaks to be drawn. Defaults to 3.
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import clugPac
#' @import GenomicRanges
#'
#' @export

patchwork_bedgraph  <- function(regions_tb,bdg_files,bdg_names,bdg_weights,epitope_colors,bdg_units="Pileup",raw_alpha=0.5,weighted_alpha=1,num_y_brks=3){
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  select  <- dplyr::select
  mutate  <- dplyr::mutate
  rename  <- dplyr::rename

  if(missing(bdg_names)){
    if(is.null(names(bdg_files))){
      names(bdg_files)  <- gsub(".bdg.gz$|.bedgraph.gz$|.bdg|.bedGraph","",basename(bdg_files))
    }
  }else{
    names(bdg_files)  <- bdg_names
  }
  if(missing(bdg_weights)){
    bdg_weights   <- rep(1,length(bdg_files))
  }else{
    if(length(bdg_weights) != length(bdg_files)){
      stop("Different number of weights/bdg files.")
    }
    if(is.null(names(bdg_weights))){
      names(bdg_weights) <- names(bdg_files)
    }
  }
  if(!missing(epitope_colors)){
    if(length(epitope_colors) != length(bdg_files)){
      epitope_colors  <- rep(epitope_colors[1],length(bdg_files))
      hide_fill_legend<- TRUE
    }
  }else{
    epitope_colors <- hue_pal()(length(bdg_files))
  }
  if(is.null(names(epitope_colors))){
    names(epitope_colors) <- names(bdg_files)
  }

  gr  <- filter(regions_tb,!trans_region) %>%
    select(seqnames1,start1,end1,strand1,region,bin_alt1_1,bin_alt1_2) %>%
    mutate(bin_alt1_1 = bin_alt1_1 - 1) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,seqnames.field = "seqnames1",start.field="start1",end.field="end1",strand.field="strand1")
  names(gr) <- gr$region

  tb_bdg  <- scan_bdg(bgz_files = bdg_files,name_vec = names(bdg_files),gr_regions = gr) %>%
    inner_join(as_tibble(gr) %>% select(region,start,end,bin_alt1_1,bin_alt1_2,strand) %>% rename(reg_start=start,reg_end=end),by="region") %>%
    mutate(width = reg_end - reg_start,
           start_frac = (start - reg_start)/width,
           end_frac = (end - reg_start)/width) %>%
    mutate_at(c("start_frac","end_frac"),
              function(x)
                case_when(x < 0 ~ 0,
                          x > 1 ~ 1,
                          TRUE ~ x)) %>%
    mutate(width=bin_alt1_2 - bin_alt1_1,
           start_bin = ifelse(strand == "+",bin_alt1_1 + width * start_frac,bin_alt1_2 - width * start_frac),
           end_bin = ifelse(strand == "+",bin_alt1_1 + width * end_frac,bin_alt1_2 - width * end_frac),
           bin_alt1_1 = ifelse(start_bin < end_bin,start_bin,end_bin),
           bin_alt1_2 = ifelse(start_bin < end_bin,end_bin,start_bin),
           weighted = score * bdg_weights[name]) %>%
    select(seqnames,start,end,strand,region,name,bin_alt1_1,bin_alt1_2,score,weighted)

  x_rng<- c(0,max(tb_bdg$bin_alt1_2))
  #y_rng<- c(0,max(tb_bdg$score))
  x_ax <- patchwork_xaxis(regions_tb)

  p <- ggplot(tb_bdg,aes(xmin=bin_alt1_1,xmax=bin_alt1_2,ymin=0,ymax=score,fill=name)) +
    facet_wrap(.~name,ncol=1,scales="free_y",strip.position = "right") +
    scale_x_continuous(name=x_ax$title,limits=x_rng,breaks=x_ax$breaks,labels=names(x_ax$breaks),expand=c(0,0)) +
    scale_y_continuous(name=bdg_units,expand=expansion(mult=c(0,0.1)),n.breaks = num_y_brks) +
    scale_fill_manual(name="Epitope",values=epitope_colors) +
    geom_rect(color=NA,alpha=raw_alpha) +
    geom_rect(color=NA,mapping=aes(ymax=weighted)) +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,-0.1,0),"lines"),
          panel.background = element_rect(fill=NA,color="black",linewidth=0.5),
          panel.grid.major.x = element_line(color="lightgray",linetype="dotted",linewidth=0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.y = unit(0,"lines"),
          strip.background = element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")
  return(p)
}
