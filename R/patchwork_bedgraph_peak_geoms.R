#' @title Patchwork bedGraph peak geoms.
#'
#' @description Returns peak annotations to be applied to patchwork_bedgraph()
#' plot which respect segment orientations etc. If an "epitope" column is
#' included, this value will respect the facets and color schemes of the same
#' column in the patchwork_bedgraph() output plot's aesthetic.
#'
#' @param regions_tb Region table, generally the output of patchwork_bins(boundaries_only=TRUE).
#' @param peak_gr GRanges object with peaks.
#' @param epitope_column Metadata column of peak_gr which differentiates peaks from different epitopes. Typically corresponds to epitope names used in patchwork_bedgraph().
#' @param respect_facets Should facetting in patchwork_bedgraph() be respected? Defaults to TRUE, otherwise peaks are drawn over all facets.
#' @param peak_alpha Transparency of each peak, defaults to 0.25.
#' @param peak_color_override Color to apply to all peaks in the even the scale_fill_manual() values in the patchwork_bedgraph() plot are to be ignored.
#' @param outline_color Color of lines to draw around each peak. Defaults to NA, i.e. outlines are not drawn.
#' @param outline_width Linewidth of outlines, if drawn. Defaults to 0.5.
#' @param outline_linetype Linetype of outlines, if drawn. Defaults to "dashed".
#'
#' @import dplyr
#' @import clugPac
#' @import tidyr
#' @import GenomicRanges
#' @import ggplot2
#'
#' @export

patchwork_bedgraph_peak_geoms <- function(regions_tb,peak_gr,epitope_column,respect_facets = TRUE,peak_alpha=0.25,peak_color_override,outline_color=NA,outline_linewidth=0.5,outline_linetype="dashed"){
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  select  <- dplyr::select
  mutate  <- dplyr::mutate
  rename  <- dplyr::rename

  gr <- peak_gr
  mcols(gr) <- NULL
  if(!missing(epitope_column)){
    if(!epitope_column %in% colnames(mcols(peak_gr))) stop("Epitope column provided, but not found in peak_gr.")
    mcols(gr)$name <- mcols(peak_gr)[,epitope_column]
  }else{
    mcols(gr)$name <- NA
  }
  ep_vals   <- unique(gr$name)

  tb_pks <- patchwork_annotations(granges_in = gr,regions_tb = regions_tb)

  if(missing(peak_color_override)){
    override_color      <- FALSE
    peak_color_override <- NA
  }
  if(nrow(tb_pks) == 0){
    g_out <- NULL
  }else{
    if(override_color & respect_facets){
      g_out <- geom_rect(data=tb_pks,mapping=aes(xmin=bin_start,xmax=bin_end,ymin=-Inf,ymax=Inf),
                         inherit.aes=FALSE,alpha=peak_alpha,fill=peak_color_override,
                         color=outline_color,linewidth=outline_linewidth,linetype=outline_linetype)
    }else if(!override_color & respect_facets){
      g_out <- geom_rect(data=tb_pks,mapping=aes(xmin=bin_start,xmax=bin_end,ymin=-Inf,ymax=Inf,fill=name),
                         inherit.aes=FALSE,alpha=peak_alpha,
                         color=outline_color,linewidth=outline_linewidth,linetype=outline_linetype)
    }else if(override_color & !respect_facets){
      tb_pks<- rename(tb_pks,epi=name)
      g_out <- geom_rect(data=tb_pks,mapping=aes(xmin=bin_start,xmax=bin_end,ymin=-Inf,ymax=Inf,fill=epi),
                         inherit.aes=FALSE,alpha=peak_alpha,
                         color=outline_color,linewidth=outline_linewidth,linetype=outline_linetype)
    }else{
      tb_pks$name <- NULL
      g_out <- geom_rect(data=tb_pks,mapping=aes(xmin=bin_start,xmax=bin_end,ymin=-Inf,ymax=Inf),
                         inherit.aes=FALSE,alpha=peak_alpha,fill=peak_color_override,
                         color=outline_color,linewidth=outline_linewidth,linetype=outline_linetype)
    }
  }
  return(g_out)
}
