#' @title
#' Read Cooler HDF5
#'
#' @description
#' Read .cool files as HDF5 using the R package rhdf5. If multiple files are
#' provided, function is run recursively on each and the results are
#' bound into a single table with a "cooler" column. Does not load entire cool
#' file into memory, but instead uses file's index to access specific rows in
#' the pixel matrices. Includes an optional caching function for larger data
#' tables by providing a TSV filename for <cache_table>. Note that integer
#' values from Cooler may be int64, which is not natively supported by R. The
#' rhdf5 package allows for the use of integer64 type using the bit64 package,
#' but as this value should only apply to genomic coordinates counts and bin ID
#' values are just converted to normal integers (there should not be 2 trillion
#' counts or bins in any normal dataset). For use with data focused on the
#' diagonal ("pyramid" plots rotated 45 degrees), the <diag_distance> argument
#' can be supplied such that pixels must be within this distance of the
#' diagonal. In addition, the X axis grange (gr_range1) is resized to be 2x the
#' diagonal distance wider, fixed about the center. This should include all
#' pixels on the left and right edges which have gr_range2 values within range,
#' but not gr_range1. The diagonal distance filter is applied AFTER fetching
#' pixels, so it is not an effective way to increase retrieval speeds. It is
#' however reflected in cached TSV files, so it will reduce file sizes.
#'
#' @param file_cool Name of .cool file(s) with the appropriate bin size. If multiple files are provided, functions is run recursively and tables are combined (including cache TSV).
#' @param gr_range1 GRange of contiguous region for X dimension. Defaults to entire genome, and must be contiguous.
#' @param gr_range2 GRange of contiguous region for Y dimension. Defaults to values of gr_range1.
#' @param diag_distace Maximum distance from the diagonal (useful for filtering pixels from tables meant for rotated plots). Filter is applied AFTER fetching pixels, so does not generally help performance.
#' @param silent Should processing time messages be suppressed? Defaults to FALSE.
#' @param max_pixels Maximum number of pixels to retrieve without erroring out. Defaults to 6.25 million (2500x2500 grid), can also be set to Inf to disregard.
#' @param cache_table Filename of TSV to cache pixel results into. If provided and this file exists, this table will be read and returned. If not found, it will be created.
#' @param overwrite_cache Should cache file be overwritten? Defaults to FALSE.
#'
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import data.table
#' @import GenomicRanges
#' @import magrittr
#' @import tictoc
#' @import clugPac
#' @import rhdf5
#' @import Rhdf5lib
#' @import data.table
#'
#' @export

read_cooler_hdf5  <- function(file_cool,gr_range1=NULL,gr_range2=NULL,diag_distance=NULL,silent=TRUE,max_pixels=6.25e6,cache_tsv=NULL,overwrite_cache = FALSE){
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  rename  <- dplyr::rename

  if(!silent) tic()
  pixels  <- NULL
  if(!is.null(cache_tsv)){
    if(file.exists(cache_tsv) & !overwrite_cache){
      #Note: MAKE SURE fread doesn't do anything with integer64, some of the math starts to chug in R. 'integer64="numeric"' ensures
      # that large integers are turned into doubles instead.
      pixels  <- fread(cache_tsv,sep="\t",header=TRUE,integer64 = "numeric") %>%
        as_tibble %>%
        mutate_at(c("seqnames1","seqnames2"), function(x) factor(x,levels=paste0('chr',c(1:22,"X","Y"))))
      if(!silent) message("Read ",prettyNum(nrow(pixels),big.mark=",")," cached pixels from ",cache_tsv,"...")
    }
  }

  if(is.null(pixels)){
    if(length(file_cool) > 1){
      if(is.null(names(file_cool))){
        names(file_cool)  <- basename(file_cool)
      }
      pixels  <- lapply(names(file_cool),function(fl_nm) {
        pxls  <- read_cooler_hdf5(file_cool = file_cool[fl_nm],
                                  gr_range1=gr_range1,
                                  gr_range2=gr_range2,
                                  diag_distance=diag_distance,
                                  silent=silent,
                                  max_pixels=max_pixels,
                                  cache_tsv=NULL) %>%
          mutate(cooler=fl_nm)
      }) %>%
        do.call(rbind,.)
    }
  }
  if(is.null(pixels)){
    #Fetch bin data.
    gr_bins   <- read_cooler_bins_hdf5(file_cool)

    #Subset relevant bins.
    if(is.null(gr_range1)){
      xbins   <- gr_bins$bin_id
      ybins   <- gr_bins$bin_id
      diag_distance <- Inf
    }else{
      #If a diagonal distance has been specified, expand gr1 to include 1x that distance on either end, otherwise left and right corners will be omitted in rotated plots.
      if(!is.null(diag_distance)){
        gr1   <- resize(gr_range1,width = width(gr_range1) + 2 * diag_distance,fix = "center")
      }else{
        gr1   <- gr_range1
        diag_distance <- Inf
      }
      #NOTE: Add one bin to the end of the x-bins list otherwise last row will happen at the start of the last block.
      xbins<- subsetByOverlaps(gr_bins,gr1)$bin_id
      xbins<- c(xbins,max(xbins+1))
      if(is.null(gr_range2)){
        ybins <- xbins
      }else{
        ybins <- subsetByOverlaps(gr_bins,gr_range2)$bin_id
        ybins <- c(ybins,max(ybins + 1))
      }
    }

    #Determine how many pixels are about to be gathered, stop if there are too many.
    pix_num   <- length(xbins) * length(ybins)
    if(pix_num > max_pixels){
      stop("Cooler subset will return up to ",prettyNum(pix_num,big.mark = ",")," pixels...to override this warning and run anyway, increase <max_pixels> or set to 'Inf'.")
    }else{
      if(!silent) message("Fetching up to ",prettyNum(pix_num,big.mark = ",")," pixels...")
    }

    #Get all bin combinations, and flip those where y < x (bottom triangle) as they are not stored in the .cool matrices.
    bin_tb  <- lapply(xbins,function(x){
      tibble(binx = x,biny=ybins)
    }) %>%
      do.call(rbind,.) %>%
      mutate(flip_x = ifelse(biny < binx,biny,NA),
             flip_y = ifelse(biny < binx,binx,NA)) %>%
      mutate(binx = ifelse(is.na(flip_x),binx,flip_x),
             biny = ifelse(is.na(flip_y),biny,flip_y))
    xbins   <- select(bin_tb,binx) %>% unlist %>% unique
    ybins   <- select(bin_tb,biny) %>% unlist %>% unique

    #Open persistent connection to HDF5.
    hdf5  <- H5Fopen(file_cool)

    #Load index.
    idx   <- h5read(hdf5,"indexes/bin1_offset")
    idx   <- tibble(bin_offset = idx) %>%
      mutate(bin_offset = as.integer(bin_offset + 1),
             bin_id = as.integer(row_number() - 1))

    #Determine which rows to retrieve.
    row_range<- filter(idx,bin_id %in% xbins) %>%
      select(bin_offset) %>%
      range

    #Fetch relevant rows in x-dimension then filter y-dimension coordinates, close connection.
    pixels  <- h5read(hdf5,name = "pixels",start=row_range[1],count=diff(row_range),bit64conversion="int")
    H5Fclose(hdf5)

    pixels  <- pixels %>%
      lapply(as.double) %>%
      as_tibble %>%
      #Filter out anything not in range 1 or range 2, including "extra" bins added for filtering step.
      filter(bin1_id %in% rev(xbins)[-1],
             bin2_id %in% rev(ybins)[-1]) %>%
      inner_join(by = c("bin1_id"="bin_id1"),
        as_tibble(gr_bins) %>% select(-width,-strand) %>% rename_all(paste0,"1")
      ) %>%
      inner_join(by = c("bin2_id"="bin_id2"),
        as_tibble(gr_bins) %>% select(-width,-strand) %>% rename_all(paste0,"2")
      ) %>%
      mutate(balanced = count * weight1 * weight2,
             log10_count = log10(count + 1)) %>%
      select(seqnames1,start_adj1,end_adj1,seqnames2,start_adj2,end_adj2,
             count,log10_count,starts_with("count"),weight1,weight2,balanced,
             start1,end1,start2,end2,bin1_id,bin2_id,everything())  %>%
      filter(abs(end2 - end1) <= diag_distance)
    if("count_wgs_nrm" %in% colnames(pixels)){
      pixels  <- mutate(pixels,log10_count_wgs_nrm = log10(count_wgs_nrm + 1))
    }
  }

  #Cache results, if needed.
  if(!is.null(cache_tsv)){
    if(!file.exists(cache_tsv) | overwrite_cache){
      if(!silent) message("Caching pixels to ",cache_tsv,"...")
      fwrite(pixels,file = cache_tsv,sep = "\t",row.names=FALSE,col.names=TRUE,quote=FALSE,)
    }
  }

  if(!silent) toc()
  return(pixels)
}
