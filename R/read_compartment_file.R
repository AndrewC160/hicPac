#' @title Read compartment eigenvector file
#' 
#' @description Read a compartment eigenvector file as output by the hic.cgp
#' pipeline and return as a granges-format tibble defining merged compartments
#' throughout the genome. Eigenvectors are reported as averages in each 
#' compartment along with standard deviations. Note that while compartments are
#' reduce()'d, empty bins are removed first, so these bins will be missing and
#' adjacent compartments may be of the same character (i.e., if an A compartment
#' contains a bin with an NA value for E1, it will be returned as two A
#' compartments flanking that missing bin).
#' 
#' @param fl_nm Filename of eigenvector file to read.
#' 
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import tibble
#' @import GenomicRanges
#'
#' @export

read_compartment_file <- function(fl_nm){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  select  <- dplyr::select
  first   <- dplyr::first
  arrange <- dplyr::arrange
  
  fread(fl_nm,sep="\t",header=TRUE) %>%
    as_tibble %>%
    filter(!is.na(E1)) %>%
    mutate(comp = ifelse(E1 > 0, "A", "B")) %>%
    group_split(comp) %>%
    setNames(sapply(.,function(tb_in) tb_in$comp[1])) %>%
    lapply(function(tb_in) {
      gr    <- makeGRangesFromDataFrame(tb_in,keep.extra.columns=TRUE)
      gr$id <- subjectHits(findOverlaps(gr,GenomicRanges::reduce(gr)))
      as_tibble(gr) %>%
        group_by(seqnames,comp,id) %>%
        summarize(start = min(start),
                  end = max(end),
                  mean_E1 = mean(E1),
                  sd_E1 = sd(E1),
                  .groups="drop") %>%
        return()
    }) %>%
    do.call(rbind,.) %>%
    mutate(sample_name = str_extract(basename(fl_nm),"^[^_]+"),
           seqnames=factor(seqnames,levels=paste0("chr",c(1:22,"X","Y")))) %>%
    arrange(seqnames,start,end) %>%
    select(seqnames,start,end,sample_name,comp,mean_E1,sd_E1) %>%
    return()
}
