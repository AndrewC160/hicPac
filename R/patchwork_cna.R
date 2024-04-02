#' @title Patchwork CNA
#'
#' @description Plot raw CNA values. Returns a GGplot object.
#'
#' @param regions_tb Region table, generally the output of patchwork_bins(boundaries_only=TRUE).
#' @param bdg_file Filename of Tabix-indexed file. Assumes index file is in place.
#' @param score_fill Fill color of rectangles to draw. If not provided, a fill based on log2(CNA raw) is applied.
#' @param num_y_brks Number of Y-axis breaks to be drawn. Defaults to 3.
#' @param line_color_in Linecolor of lines separating segments. Defaults to gray20.
#' @param linewidth_in Linewidth of lines separating segments. Defaults to 0.5.
#' @param linetype_in Linetype of lines separating segments. Defaults to "dotted."
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import clugPac
#' @import ggsci
#' @import GenomicRanges
#'
#' @export

patchwork_cna <- function(regions_tb,bdg_file,score_fill,num_y_brks=3,line_color_in="gray20",linewidth_in = 0.5,linetype_in="dotted"){
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  select  <- dplyr::select
  mutate  <- dplyr::mutate
  rename  <- dplyr::rename

  gr  <- filter(regions_tb,!trans_region) %>%
    select(seqnames1,start1,end1,strand1,region,bin_alt1_1,bin_alt1_2) %>%
    mutate(bin_alt1_1 = bin_alt1_1 - 1) %>%
    arrange(bin_alt1_1,bin_alt1_2) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,seqnames.field = "seqnames1",start.field="start1",end.field="end1",strand.field="strand1")
  names(gr) <- gr$region

  tb_bdg  <- lapply(names(gr), function(nm) {
    scan_bdg(bdg_file,gr_regions = gr[nm]) %>%
      as_tibble
  }) %>% do.call(rbind,.) %>%
    full_join(as_tibble(gr) %>% select(region,start,end,bin_alt1_1,bin_alt1_2,strand) %>% rename(reg_start=start,reg_end=end),by="region") %>%
    group_by(region) %>%
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
           bin_alt1_2 = ifelse(start_bin < end_bin,end_bin,start_bin)) %>%
    select(seqnames,start,end,strand,region,name,bin_alt1_1,bin_alt1_2,score) %>%
    ungroup

  x_rng<- c(0,max(tb_bdg$bin_alt1_2))
  x_ax <- patchwork_xaxis(regions_tb)

  if(missing(score_fill)){
    rect_geoms  <- geom_rect(color=NA)
  }else{
    rect_geoms  <- geom_rect(color=NA,fill=score_fill)
  }

  p <- ggplot(tb_bdg,aes(xmin=bin_alt1_1,xmax=bin_alt1_2,ymin=0,ymax=score,fill=score)) +
    scale_x_continuous(name=x_ax$title,limits=x_rng,breaks=x_ax$breaks,labels=names(x_ax$breaks),expand=c(0,0)) +
    scale_y_continuous(name="CNA\nraw",expand=expansion(mult=c(0,0.1)),n.breaks = num_y_brks,
                       labels=function(x) ifelse(x==0,"0\n",paste0(x/1e3,"k"))) +
    scale_fill_viridis() +
    rect_geoms +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,-0.1,0),"lines"),
          panel.background = element_rect(fill=NA,color="black",linewidth=0.5),
          panel.grid.major.x = element_line(color=line_color_in,linetype=linetype_in,linewidth=linewidth_in),
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
