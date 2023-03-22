#' @title Patchwork annotations.
#'
#' @description Given a GenomicRanges object of genome annotations and a
#' region table from patchwork_bins(boundaries_only=TRUE),subset and clip the
#' contents of these ranges to identify features which fall within each contig
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

patchwork_annotations <- function(granges_in,regions_tb=NULL,hic_file=NULL,gr_segments=NULL){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange

  if(is.null(regions_tb)){
    if(is.null(hic_file) | is.null(gr_segments)) stop("Either a regions table from patchwork_bins() is required or both the hic_file and gr_segments arguments must be provided.")
    regions_tb <- patchwork_bins(gr_list = gr_segments,hic_file = fl_hic)
  }
  gr_reg  <- filter(regions_tb,!trans_region) %>%
    makeGRangesFromDataFrame(seqnames.field = "seqnames1",
                             start.field = "start1",
                             end.field = "end1",
                             strand.field = "strand1",
                             keep.extra.columns = TRUE,
                             ignore.strand = FALSE)

  lapply(c(1:length(gr_reg)), function(i) {
    gr <- gr_reg[i]
    reg<- gr$region
    rv <- as.character(strand(gr)) == "-"
    x_rng     <- c(start(gr),end(gr))
    x_rng_bin <- c(gr$bin_alt1_1,gr$bin_alt1_2)

    clugPac::clip_granges(granges_in,gr,return_as_tibble = TRUE,include_clipped_col = TRUE) %>%
      rowwise %>%
      mutate(reversed = rv,
             start_frac = (start - x_rng[1]) / diff(x_rng),
             end_frac = (end - x_rng[1]) / diff(x_rng),
             start_frac = ifelse(reversed,1-start_frac,start_frac),
             end_frac = ifelse(reversed,1-end_frac,end_frac),
             strand = as.character(strand),
             strand = case_when(reversed & strand == "+" ~ "-",
                                reversed & strand == "-" ~ "+",
                                TRUE ~ strand),
             trunc_start = ifelse(reversed,start_clipped != 0,end_clipped != 0),
             trunc_end = ifelse(reversed,end_clipped != 0,start_clipped != 0),
             bin_start = x_rng_bin[1] + start_frac * diff(x_rng_bin),
             bin_end = x_rng_bin[1] + end_frac * diff(x_rng_bin),
             region = reg,
             start_flip = ifelse(bin_end > bin_start,bin_end,bin_start),
             end_flip = ifelse(bin_end > bin_start, bin_start,bin_end),
             bin_start = start_flip,
             bin_end = end_flip) %>%
      select(-start_flip,-end_flip) %>%
      mutate_at(c("bin_start","bin_end"), function(x)
        case_when(x < x_rng_bin[1] ~ x_rng_bin[1],
                  x > x_rng_bin[2] ~ x_rng_bin[2],
                  TRUE ~ x)) %>%
      ungroup %>%
      select(-start_clipped,-end_clipped,-width,-start_frac,-end_frac) %>%
      mutate(bin_start = ifelse(bin_start < x_rng_bin[1],x_rng_bin[1],bin_start),
             bin_end = ifelse(bin_end > x_rng_bin[2],x_rng_bin[2],bin_end)) %>%
      return()
  }) %>%
    do.call(rbind,.) %>%
    return()
}
