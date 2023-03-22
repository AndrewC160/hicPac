#' @title Patchwork bins
#'
#' @description
#' Get bins and coordinates for a given contig of DNA starting with a list of
#' genomic ranges of DNA segments which comprise that contig. Can be named, but
#' if not each segment will be labeled "reg1", "reg2", etc. Will also append a
#' "bin_alt" coordinate system which is specific to this contig and is in units
#' of pixels (the start of the contig is pixel 1, aka bin_alt1_1). Direction of
#' each segment is indicated by its range's strand ("-" indicates 3' to 5'
#' direction).
#'
#' @param gr_list List of GRanges to combine.
#' @param hic_file Cooler filename from which to retrieve bins.
#' @param boundaries_only Defaults to TRUE, in which case only bins which bound each segment are retrived. If FALSE, all bins are retrieved.
#'
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import magrittr
#'
#' @export

patchwork_bins    <- function(gr_list,hic_file,boundaries_only=TRUE){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange

  if(inherits(gr_list,"GRanges")){
    gr_list$region <- paste0("reg",c(1:length(gr_list)))
    gr_list <- split(gr_list,gr_list$region)
  }
  tb_bins <- read_cooler_bins_hdf5(file_cooler = hic_file,granges_list = gr_list) %>% as_tibble
  bin_ids <- lapply(tb_bins$bin_alt,function(x){
    tibble(binx = x,biny=tb_bins$bin_alt)
  }) %>%
    do.call(rbind,.) %>%
    filter(binx <= biny) #Only take upper triangle.

  tb_out  <-
    cbind(as_tibble(tb_bins[bin_ids$binx,]) %>% select(-width) %>% rename_all(paste0,"1"),
          as_tibble(tb_bins[bin_ids$biny,]) %>% select(-width) %>% rename_all(paste0,"2")) %>%
    as_tibble %>%
    mutate(trans_region = region1 != region2,
           region = ifelse(!trans_region,region1,paste0(region1,"-",region2))) %>%
    select(-region1,-region2) %>%
    #Reversed coordinates will flip x/y bins (i.e. reflect to the bottom triangle), so switch them back.
    mutate(flip_x = ifelse(bin_id1 > bin_id2,bin_id2,bin_id1),
           flip_y = ifelse(bin_id1 > bin_id2,bin_id1,bin_id2),
           bin_id1 = flip_x,
           bin_id2 = flip_y) %>%
    select(-flip_x,-flip_y)
  if(boundaries_only){
    tb_out <- tb_out %>%
      group_by(region,trans_region) %>%
      summarize(seqnames1 = first(seqnames1),
                strand1 = first(strand1),
                start1 = min(start1),
                end1 = max(end1),
                bin_id1_1 = min(bin_id1),
                bin_id1_2 = max(bin_id1),
                bin_alt1_1= min(bin_alt1),
                bin_alt1_2= max(bin_alt1),
                seqnames2 = first(seqnames2),
                strand2 = first(strand2),
                start2 = min(start2),
                end2 = max(end2),
                bin_id2_1 = min(bin_id2),
                bin_id2_2 = max(bin_id2),
                bin_alt2_1= min(bin_alt2),
                bin_alt2_2= max(bin_alt2),
                .groups = "drop")
  }
  return(tb_out)
}
