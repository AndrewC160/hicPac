#' @title
#' Dump cooler
#'
#' @description
#' Dump cooler contents from a given region, either in one square (same range
#' in two dimensions) or in an asymmetric rectangle (one range for X dimension,
#' another for Y).
#'
#' @param gr_range1 GRange of contiguous region for X dimension. Defaults to entire genome.
#' @param gr_range2 GRange of contiguous region for Y dimension. Defaults to values of gr_range1.
#' @param file_cool Filename of .cool file with the appropriate bin size.
#' @param exe_cooler Filename of cooler executable. Defaults to '~/anaconda3/envs/EagleC/bin/cooler'.
#'
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import data.table
#' @import GenomicRanges
#' @import magrittr
#' @import scales
#' @import ggplot2
#' @import tictoc
#' @import clugPac
#'
#' @export

dump_cooler   <- function(gr_range1=NULL,gr_range2=NULL,file_cool,silent=FALSE,exe_cooler="/N/u/aclugston/anaconda3/envs/EagleC/bin/cooler"){
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  if(is.null(gr_range1)){
    gr_1  <- ""
    gr_1_desc <- ""
    gr_2  <- ""
    gr_2_desc <- ""
  }else{
    gr_1  <- paste("--range",as.character(gr_range1))
    gr_1_desc <- grange_desc(gr_range1)
    if(is.null(gr_range2)){
      gr_2<- paste("--range2",as.character(gr_range1))
      gr_2_desc <- grange_desc(gr_range1)
    }else{
      gr_2<- paste("--range2",as.character(gr_range2))
      gr_2_desc <- grange_desc(gr_range2)
    }
  }
  seq_adjs<- get_seqsizes_adj()
  if(!silent) tic()
  pixels<- paste(exe_cooler,"dump --table pixels",gr_1,gr_2,
                 "--header --one-based-ids --one-based-starts --join --annotate weight --balanced",
                 file_cool) %>%
    fread(cmd=.,sep="\t",header = TRUE) %>%
    as_tibble %>%
    mutate(seqnames1=factor(chrom1,levels=paste0("chr",c(1:22,"X","Y"))),
           seqnames2=factor(chrom2,levels=paste0("chr",c(1:22,"X","Y"))),
           adj1 = seq_adjs[as.character(seqnames1)],
           adj2 = seq_adjs[as.character(seqnames2)],
           start1_adj = start1 + adj1,
           end1_adj = end1 + adj1,
           start2_adj = start2 + adj2,
           end2_adj = end2 + adj2) %>%
    select(seqnames1,start1_adj,end1_adj,seqnames2,start2_adj,end2_adj,count,weight1,weight2,balanced,start1,end1,start2,end2) %>%
    group_by(start1) %>%
    mutate(bin1=cur_group_id()) %>%
    ungroup %>%
    group_by(start2) %>%
    mutate(bin2=cur_group_id()) %>%
    ungroup %>%
    mutate(log10_count = log10(count))
  if(!silent) toc()
  return(pixels)
}
