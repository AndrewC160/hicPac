#' @title
#' Patchwork matrix.
#'
#' @description
#' Retrive bins and pixels of HiC data corresponding to a contig assembled from
#' a series of genomic regions. Starts with a list of genomic regions (gr_list)
#' from 5' to 3' whose direction is indicated by strand, and a hic Cooler file
#' (hic_file) from which to retrieve data.
#'
#' Given the instability of some rhdf5 functions (currently), this function
#' requires the use of a cache system. First, a temporary directory is saved as
#' a file with ".tsv" replaced with "_tmp" and every region is read directly
#' into a TSV file in this directory. Once all TSVs have been found, these are
#' read using data.table::fread() and concatenated into the output table. Once
#' the table's existence is confirmed, the temporary directory and its contents
#' are recursively unlinked. Hopefully this will turn out to be excessive once
#' I figure out the cause of so many rhdf5 failures, but for now this will do.
#' Worst case scenario, files can be "ratcheted through" such that following
#' each crash, the process can be resumed at the last incomplete file.
#'
#' @param gr_list List of GRanges (can be named) which will be retrieved, with ranges coming in order of 5' to 3' in the contig in question. Strand indicates a given segment's direction.
#' @param hic_file Cooler filename from which to retrieve bins and pixel data.
#' @param pixel_tsv TSV file in which to cache output. If <over_write> is FALSE and <pixel_tsv> exists, this file will be read from disk instead of retrieving from the Cooler file.
#' @param over_write Should existing <pixel_tsv> file be overwritten? Defaults to FALSE.
#' @param max_pixels Maximum number of pixels to retrieve for a given segment, defaults to 6.25E6. Can be increased to Inf.
#'
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import GenomicRanges
#' @import magrittr
#'
#' @export

patchwork_matrix  <- function(gr_list,hic_file,pixel_tsv,over_write=FALSE,max_pixels = 6.25E6){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange

  if(missing(pixel_tsv)) stop("A cache file is required (pixel_tsv).")
  if(is.null(names(gr_list))) names(gr_list) <- paste0("reg_",c(1:length(gr_list)))
  if(missing(pixel_tsv)) pixel_tsv <- ""
  if(file.exists(pixel_tsv) & !over_write){
    tb_pix  <- fread(pixel_tsv,integer64 = "double",sep="\t") %>%
      as_tibble %>%
      mutate_at(c("seqnames1","seqnames2"),
                function(x) factor(x,levels=paste0("chr",c(1:22,"X","Y"))))
  }else{
    tb_bins <- patchwork_bins(gr_list,hic_file,boundaries_only = FALSE)
    tb_regs   <- tb_bins %>%
      group_by(region) %>%
      summarize(seqnames1 = first(seqnames1),
                strand1 = first(strand1),
                start1 = min(start1),
                end1 = max(end1),
                seqnames2 = first(seqnames2),
                strand2 = first(strand2),
                start2 = min(start2),
                end2 = max(end2),
                .groups = "drop")

    #Read pixel data and merge into bins.
    ls_pix <- setNames(vector(mode = "list",length = nrow(tb_regs)),tb_regs$region)
    #Create temporary directory for progress.
    tmp_dir<- gsub(".tsv","_tmp",pixel_tsv)
    dir.create(tmp_dir,showWarnings = FALSE)
    for(i in 1:length(ls_pix)){
      tb    <- tb_regs[i,]
      reg_nm<- tb$region
      tmp_fl<- file.path(tmp_dir,paste0(reg_nm,".tmp"))
      grA <- makeGRangesFromDataFrame(tb,
                                      seqnames.field = "seqnames1",
                                      start.field = "start1",
                                      end.field = "end1")
      grB <- makeGRangesFromDataFrame(tb,
                                      seqnames.field = "seqnames2",
                                      start.field = "start2",
                                      end.field = "end2")
      ls_pix[[i]] <- read_cooler_hdf5(file_cool = hic_file,
                                    cache_tsv = tmp_fl,
                                    gr_range1 = grA,
                                    gr_range2 = grB,
                                    max_pixels = max_pixels,
                                    return_table=FALSE)
    }
    tb_pix <- lapply(names(ls_pix), function(nm) {
      fread(ls_pix[[nm]],header=TRUE) %>%
        as_tibble %>%
        mutate(region=nm)
    }) %>%
      do.call(rbind,.) %>%
      right_join(by=c("bin1_id"="bin_id1","bin2_id"="bin_id2","region"="region"),multiple="all",
                 select(tb_bins,strand1,strand2,bin_id1,bin_id2,bin_alt1,bin_alt2,region)
      ) %>%
      drop_na %>%
      arrange(region,bin_alt1,bin_alt2) %>%
      group_by(region) %>%
      mutate(boundary = row_number() == 1 | row_number() == n()) %>%
      ungroup %>%
      filter(count > 0 | boundary) %>%
      mutate(distance = abs(bin_alt2 - bin_alt1))
    fwrite(tb_pix,file = pixel_tsv,sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)

    #If everything worked, remove the temporary directory.
    if(file.exists(pixel_tsv)) unlink(tmp_dir,recursive=TRUE)
  }
  return(tb_pix)
}
