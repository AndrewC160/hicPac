#' @title
#' Read cooler genomic bins
#'
#' @description
#' Given a cooler file, retrieve genomic bin coordinates and return as a
#' GRanges object. Does not open file continually, but requires that there not
#' be an existing connection to another HDF5 based on rhdf5 requirements. If
#' strand of input GRanges is negative, this is assumed to mean that said
#' region is in the reverse direction, and so the bin's strand value will be
#' matched.
#'
#' @param file_cooler Cooler file from which to retrive bin information.
#' @param granges_list Genomic ranges to be incorporated. If one or more is provided, bins will be retrieved for each, plus an additional arbitrary bin ID.
#' @param ignore_strand Should strand orientation be ignored among input granges?
#' @param disable_file_lock Should the locking of HDF5 files be disabled? Defaults to FALSE, but can be set to TRUE if multiple instances will be accessing the same Cooler simultaneously.
#'
#'
#' @import rhdf5
#' @import dplyr
#' @import tibble
#' @import magrittr
#' @import GenomicRanges
#' @import clugPac
#'
#' @export

read_cooler_bins_hdf5 <- function(file_cooler,granges_list = NULL,ignore_strand=FALSE,disable_file_lock=FALSE){
  sq_offs <- get_seqsizes_adj()
  if(disable_file_lock){
    h5disableFileLocking()
  }else{
    h5enableFileLocking()
  }

  hdf5    <- H5Fopen(file_cooler,flags="H5F_ACC_RDONLY")
  bin_data<- h5read(hdf5,name="bins",bit64conversion="double")
  H5Fclose(hdf5)
  gr_bins <- bin_data %>%
    as_tibble %>%
    mutate(bin_id = as.integer(row_number()-1),
           start_adj = start + sq_offs[as.character(chrom)],
           end_adj = end + sq_offs[as.character(chrom)]) %>%
    mutate(chrom = factor(chrom,levels=paste0("chr",c(1:22,"X","Y")))) %>%
    filter(!is.na(chrom)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  if(!is.null(granges_list)){
    if(!inherits(granges_list,"GRangesList")){
      granges_list <- GRangesList(granges_list)
    }
    if(is.null(names(granges_list))){
      names(granges_list) <- paste0("reg",c(1:length(granges_list)))
    }
    gr_bins <- lapply(names(granges_list), function(nm){
      gr  <- granges_list[[nm]]
      rvrs<- as.character(strand(gr)) == "-" | ignore_strand
      sort(subsetByOverlaps(gr_bins,gr),decreasing = rvrs) %>%
        as_tibble %>%
        mutate(region = nm,
               strand = ifelse(rvrs,"-","+"))
    }) %>%
      do.call(rbind,.) %>%
      mutate(bin_alt = row_number()) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                               seqnames.field = "seqnames",
                               start.field = "start",
                               end.field = "end",
                               strand.field = "strand")
  }
  return(gr_bins)
}
