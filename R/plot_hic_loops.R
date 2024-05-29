#' @title Plot HiC loops
#'
#' @description Given a pixel tibble (non-rotated) with at least one column with
#' score values, draw geom_curve()s between pixel bins with color and/or alpha
#' values tied to that score column. If 'trap' regions are provided via gr_trap,
#' illustrated loops will be restricted to those which start/end in one of those
#' trap windows. For instance, a 10kb trap centered on the promoter region for a
#' gene will restrict loops to those which start/end on the gene promoter region
#' and thus may indicate promoter/enhancer loops.
#'
#' Note that loops are only drawn if both ends are within the x-axis range; no
#' loops to outside regions are shown. If x_rng is specified, loops to bins off
#' of the x-range will be dropped. Loops cannot be drawn between the same bin,
#' so pixels on the diagonal will not be drawn (distance must be greater than
#' zero).
#'
#' Also note that color_col_name and alpha_col_name default to NULL, in which
#' case color and alpha scales (respectively) will be drawn without a legend,
#' but the scales will be unaffected.
#'
#' @param pixel_tibble Tibble containing pixel data (not rotated) with the score values of interest.
#' @param gr_trap GRanges object of trap regions to filter by. Defaults to NULL, in which case all loops are shown.
#' @param gr_window GRanges object denoting the x-range window to be drawn. Defaults to NULL, in which case the window will be bounded by the range of X-values in the data.
#' @param invert_y Should the Y-axis be inverted (loops open upward)? Defaults to TRUE.
#' @param color_col Column name to which color values should be tied. Defaults to "log10_counts".
#' @param color_col_name Name of color column to be used on the legend. Defaults to "-log10(Count)".
#' @param alpha_col Column name to which alpha values should be tied. Defaults to NULL.
#' @param alpha_col_name Name of alpha column to be used on the legend. Defaults to NULL.
#' @param default_color Default color in the event no color scale is used. Defaults to "black".
#' @param default_alpha Default alpha values in the event no alpha scale is used. Defaults to 1.
#'
#' @import ggplot2
#' @import tibble
#' @import GenomicRanges
#' @import magrittr
#'
#' @export

plot_hic_loops <- function(pixel_tibble,gr_trap=NULL,gr_window=NULL,invert_y=TRUE,
                           color_col=NULL,color_col_name=NULL,
                           alpha_col=NULL,alpha_col_name=NULL,
                           default_color="black",default_alpha=1){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  start   <- GenomicRanges::start
  end     <- GenomicRanges::end
  between <- data.table::between

  tb_p <- pixel_tibble %>%
    mutate(d_col = abs(start1 - start2),
           w_col = end1-start1) %>%
    filter(d_col > w_col,
           d_col > 0) %>%
    mutate(x=(start1 + end1)/2,
           xend=(start2 + end2)/2,
           nlog10Q = -log10(qvalue)) %>%
    rename(seqnames = seqnames1) %>%
    select(-seqnames2,-starts_with("start"),-starts_with("end"),-d_col,-w_col)

  color_legend_mode <- "none"
  if(!is.null(color_col)){
    tb_p <- mutate(tb_p,score_color := !!as.name(color_col))
    scale_color <- scale_color_gradient(name=color_col_name,low="blue",high="red",na.value = "green")
    if(!is.null(color_col_name)){
      color_legend_mode <- NULL
    }
  }else{
    tb_p <- mutate(tb_p,score_color=1)
    scale_color <- scale_color_gradient(low=default_color,high=default_color)
  }

  alpha_legend_mode <- "none"
  if(!is.null(alpha_col)){
    tb_p <- mutate(tb_p,score_alpha := !!as.name(alpha_col))
    scale_alpha <- scale_alpha_continuous(name=alpha_col_name)
    if(!is.null(alpha_col_name)) alpha_legend_mode <- NULL
  }else{
    tb_p <- mutate(tb_p,score_alpha=1)
    scale_alpha <- scale_alpha_continuous(range=c(default_alpha,default_alpha))
  }

  if(is.null(gr_window)){
    x_rng <- c(min(tb_p$x),max(tb_p$xend))
  }else{
    x_rng <- c(start(gr_window),end(gr_window))
  }
  gr_wnd<- GRanges(as.character(tb_p$seqnames[1]),IRanges(x_rng[1],x_rng[2]))
  tb_p <- filter(tb_p,
            between(x,start(gr_wnd),end(gr_wnd),incbounds = FALSE) &
            between(xend,start(gr_wnd),end(gr_wnd),incbounds = FALSE))

                 # between(x,x_rng[1],x_rng[2],incbounds = FALSE) &
                 #      between(xend,x_rng[1],x_rng[2],incbounds = FALSE))

  if(!is.null(gr_trap)){
    c_keep  <- rep(FALSE,nrow(tb_p))
    geom_trap_list <- vector(mode = "list",length = length(gr_trap))
    for(i in 1:length(gr_trap)){
      gr_t <- gr_trap[i]
      t_keep <- between(tb_p$x,start(gr_t),end(gr_t),incbounds = TRUE) | between(tb_p$xend,start(gr_t),end(gr_t),incbounds = TRUE)
      c_keep <- c_keep | t_keep
      geom_trap_list[[i]]  <- geom_vline(xintercept = c(start(gr_t),end(gr_t)),color="red")
    }
    tb_p <- tb_p[c_keep,]
  }else{
    geom_trap_list <- NULL
  }

  if(nrow(tb_p) == 0){
    p <- ggplot() + geom_blank() +
      scale_x_continuous(name=grange_desc(gr_wnd),expand=c(0,0),limits=x_rng,
                         labels = prettyBP) +
      scale_y_reverse(limits=c(1,0)) +
      geom_trap_list
  }else{
    if(invert_y){
      p <- ggplot(tb_p,aes(x=x,xend=xend,y=0,yend=0,color=score_color,alpha=score_alpha)) +
        scale_y_reverse(expand=c(0,0),limits=c(1,0))
    }else{
      p <- ggplot(tb_p,aes(x=xend,xend=x,y=0,yend=0,color=score_color,alpha=score_alpha)) +
        scale_y_reverse(expand=c(0,0),limits=c(0,-1))
    }

    p <- p +
      facet_wrap(.~cooler,ncol=1,strip.position = "right",scales="fixed") +
      scale_x_continuous(name=grange_desc(gr_wnd),expand=c(0,0),limits=x_rng,labels = prettyBP) +
      scale_color +
      scale_alpha +
      geom_curve() +
      geom_trap_list +
      guides(color=color_legend_mode,alpha=alpha_legend_mode)
  }
  p <- p +
    theme(plot.background = element_rect(fill="white",color=NA),
          plot.margin = unit(c(0,0,0,0),"lines"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill=NA,color="black",linewidth=0.5),
          panel.spacing.y = unit(0.05,"lines"),
          strip.background = element_blank(),
          strip.text.y.right = element_text(angle=0),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  return(p)
}
