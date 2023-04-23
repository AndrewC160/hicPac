#' @title
#' Get axis labels
#'
#' @description
#' Given a tibble returned by read_cooler_hdf5(), determine labels for x/y axes
#' and return as a table or vector. Specify <break_num> to approximate how many
#' breaks are required, and all chromosome starts will be included. If more
#' chromosomes are contained in the data than there are breaks, only chromosome
#' start points will be annotated, but they all will be regardless of break_num.
#'
#' Typical format: "Chr1 | 1,234,567 | 106,434,567 | Chr2 | 1,234,567"
#'
#' @param pix_in Matrix tibble in the format returned by read_matrix_file().
#' @param axis_in Which axis to return labels for. Defaults to "x", can also be "y". For distance in a rotated matrix, use "y_hic".
#' @param break_num Number of breaks to approximate. All chromosome starts must be annotated, however, so if there are more chromosomes than <break_num> it will be effectively ignored.
#' @param out_as_vector Should the output table be simplified down to a named vector with the offset start positions as values? Defaults to true.
#'
#' @import dplyr
#' @import magrittr
#' @import GenomicRanges
#'
#' @export

get_axis_labs <- function(pix_in,axis_in = "x",break_num = 5,out_as_vector=TRUE){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  if(axis_in == "y_hic"){
    if(!"y_coords" %in% colnames(pix_in)){
      stop("No 'y_coords' column found in the pixel table, does it still need to be rotated?")
    }else{
      y_rng   <- 2 * (range(pix_in$y_coords) - min(pix_in$y_coords))
      y_seq   <- c(seq(y_rng[1],y_rng[2],length.out=break_num + 1)[-1])
      y_brks  <- (y_seq)/2 + min(pix_in$y_coords)
      names(y_brks) <- prettyBP(y_seq)
    }
    brks_out  <- y_brks
  }else{
    if(axis_in == "y"){
      mtx_out  <- pix_in %>%
        rename(start = start2,
               start_adj = start_adj2,
               seqnames = seqnames2)
    }else{
      mtx_out  <- pix_in %>%
        rename(start = start1,
               start_adj = start_adj1,
               seqnames = seqnames1)
    }
    mtx_out   <- select(mtx_out,seqnames,start,start_adj)
    seq_num   <- length(unique(mtx_out$seqnames))

    if(seq_num > break_num){
      brks_out<- mtx_out %>%
        group_by(seqnames) %>%
        summarize(start = 1,
                  start_adj = min(start_adj),
                  .groups = "drop") %>%
        rename(label = seqnames)
    }else if(seq_num == 1){
      seq_offs<- get_seqsizes_adj()
      brks    <- seq(min(mtx_out$start),
                     max(mtx_out$start) - (max(mtx_out$start) - min(mtx_out$start)) / break_num,
                     length.out=break_num)

      brks_out<- mtx_out %>%
        mutate(break_num = cut(start,breaks = brks,labels=FALSE,include.lowest = TRUE)) %>%
        group_by(break_num,seqnames) %>%
        summarize(start = min(start),
                  start_adj = min(start_adj),
                  .groups="drop") %>%
        select(-break_num) %>%
        mutate(seq_label = ifelse(as.character(seqnames) == lag(as.character(seqnames),default = ""),
                                  "",prettyTitle(as.character(seqnames))),
               pos_label = ifelse(start == 0,"",prettyNum(start,big.mark=",")),
               label = case_when(pos_label == "" & seq_label != "" ~ seq_label,
                                 pos_label != "" & seq_label == "" ~ pos_label,
                                 TRUE ~ paste0(seq_label,":",pos_label))) %>%
        select(-seq_label,-pos_label)
    }
    if(out_as_vector){
      brks_out  <- vectify(brks_out,start_adj,label)
    }
  }
  return(brks_out)
}
