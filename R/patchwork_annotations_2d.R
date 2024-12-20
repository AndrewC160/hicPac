#' @title Patchwork annotations 2D.
#'
#' @description Given a GenomicRanges object of genome annotations representing 
#' two coordinates in a patchwork (for example, a translocation), and a region 
#' table from patchwork_bins(boundaries_only=TRUE), subset and clip the contents
#' of these ranges to identify features which fall within the matrix.
#' 
#' 
#' segment's range. Annotates truncated features with "trunc_start" and
#' "trunc_end" columns, and while start/end coordinates are provided, X and Y
#' coordinates for plotting are found in "bin_start" and "bin_end" columns,
#' which have contig-specific values. These values also respect reversal:
#' features which fall on segments that have been added in the 3'-5' direction
#' are reversed ("start" and "end" coordinates run backwards), and strand
#' information is inverted ("+" strand genes are labeld "-"). If <regions_tb>
#' is not provided, <hic_file> and <gr_segments> can be substituted to execute
#' patchwork_bins() within the function.
#'
#' @param granges_in GenomicRanges object containing annotations to transform.
#' @param regions_tb Table from patchwork_bins(boundaries_only=TRUE) defining the contig in question.
#' @param hic_file Cooler file from which to extract bin information if <regions_tb> is NULL; ignored otherwise.
#' @param gr_segments GenomicRanges object containing contig segment information to be used if <regions_tb> is NULL; ignored otherwise.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import GenomicRanges
#' @import clugPac
#'
#' @export

# tb_annots <- patchwork_annotations(makeGRangesFromDataFrame(tb_svs,keep.extra.columns = TRUE,seqnames.field = "seqnames1",start.field="pos1",end.field="pos1"),
#                                    regions_tb = tb_regs)
# granges1_in <- makeGRangesFromDataFrame(tb_svs,keep.extra.columns = FALSE,seqnames.field = "seqnames1",start.field="pos1",end.field="pos1")
# granges2_in <- makeGRangesFromDataFrame(tb_svs,keep.extra.columns = FALSE,seqnames.field = "seqnames2",start.field="pos2",end.field="pos2")

patchwork_annotations_2d <- function(granges1_in,granges2_in,regions_tb,hic_file,gr_segments){
  if(length(granges1_in) != length(granges2_in)){
    stop("2D GRanges must be of the same length; currently granges 1 and 2 are ",paste0(length(granges1_in)," and ",length(granges2_in),", respectively."))
  }
  if("index_column" %in% c(colnames(as_tibble(granges1_in)),
                           colnames(as_tibble(granges2_in)))) stop("'index_column' cannot be included as a column name for either input GRanges 1 or 2.")
  gr1 <- granges1_in
  gr2 <- granges2_in
  gr1$index_column <- paste0("idx_",1:length(gr1))
  gr2$index_column <- paste0("idx_",1:length(gr2))
  
  anno1 <- patchwork_annotations(gr1,regions_tb=tb_regs)
  anno2 <- patchwork_annotations(gr2,regions_tb=tb_regs)
  
  inner_join(anno1,anno2,by="index_column",suffix = c("1","2")) %>%
    select(-index_column) %>%
    return()
}

