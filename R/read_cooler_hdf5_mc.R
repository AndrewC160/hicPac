#' @title
#' Read Cooler HDF5 (Multi-core)
#'
#' @description
#' Read .cool files as HDF5 using the R package rhdf5, and break up the job into multiple cores.If
#' multiple files are provided, function is run recursively on each and the results are bound into a
#' single table with a "cooler" column. Does not load entire cool file into memory, but instead uses
#' file's index to access specific rows in the pixel matrices. Includes an optional caching function
#' for larger data tables by providing a TSV filename for <cache_table>. Note that integer values
#' from Cooler may be int64, which is not natively supported by R. The rhdf5 package allows for the
#' use of integer64 type using the bit64 package, but as this value should only apply to genomic
#' coordinates counts and bin ID values are just converted to normal integers (there should not be
#' 2 trillion counts or bins in any normal dataset). For use with data focused on the diagonal
#' ("pyramid" plots rotated 45 degrees), the <diag_distance> argument can be supplied such that
#' pixels must be within this distance of the diagonal. In addition, the X axis grange (gr_range1)
#' is resized to be 2x the diagonal distance wider, fixed about the center. This should include all
#' pixels on the left and right edges which have gr_range2 values within range, but not gr_range1.
#' The diagonal distance filter is applied AFTER fetching pixels, so it is not an effective way to
#' increase retrieval speeds. It is however reflected in cached TSV files, so it will reduce file
#' sizes.
#'
#' In some cases, it is preferable to produce the TSV cache file and NOT load its contents into
#' memory. In these cases, <return_table> can be set to FALSE and the pixel TSV file will be
#' returned; in the event this file exists it will not be read into memory.
#'
#' Note that for multiple tasks to be completed at once, file locking must be disabled and thus
#' this function is run with the equivalent of read_cooler_hdf5(disable_file_lock=TRUE).
#'
#' @param file_cool Name of .cool file(s) with the appropriate bin size. If multiple files are provided, functions is run recursively and tables are combined (including cache TSV).
#' @param gr_range1 GRange of contiguous region for X dimension. Defaults to entire genome, and must be contiguous.
#' @param gr_range2 GRange of regions for Y dimension. Defaults to values of gr_range1, and can be multiple ranges.
#' @param diag_distace Maximum distance from the diagonal (useful for filtering pixels from tables meant for rotated plots). Filter is applied AFTER fetching pixels, so does not generally help performance.
#' @param silent Should processing time messages be suppressed? Defaults to FALSE.
#' @param max_pixels Maximum number of pixels to retrieve without erroring out. Defaults to 6.25 million (2500x2500 grid), can also be set to Inf to disregard.
#' @param cache_table Filename of TSV to cache pixel results into. If provided and this file exists, this table will be read and returned. If not found, it will be created.
#' @param overwrite_cache Should cache file be overwritten? Defaults to FALSE.
#' @param return_table Should the table be returned? Defaults to TRUE, but if the goal is to produce the TSV cache file this option can avoid extra reading and memory requirements. Requires a cache_tsv file, if FALSE.
#' @param n_cores Number of cores to process data with; defaults to 10.
#' @param n_chunks Number of chunks to process data in; defaults to 20.
#'
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import data.table
#' @import GenomicRanges
#' @import IRanges
#' @import magrittr
#' @import tictoc
#' @import clugPac
#' @import rhdf5
#' @import Rhdf5lib
#' @import data.table
#' @import parallel
#'
#' @export

read_cooler_hdf5_mc  <- function(file_cool,gr_range1=NULL,gr_range2=NULL,diag_distance=NULL,silent=TRUE,max_pixels=6.25e6,cache_tsv=NULL,overwrite_cache = FALSE,return_table=TRUE,n_cores=10,n_chunks=20){
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  filter  <- dplyr::filter
  rename  <- dplyr::rename

  pixels  <- NULL
  if(!is.null(cache_tsv)){
    if(file.exists(cache_tsv) & !overwrite_cache){
      # Note: MAKE SURE fread doesn't do anything with integer64, some of the math starts to chug in R. 'integer64="numeric"' ensures
      # that large integers are turned into doubles instead.
      if(return_table){
        pixels  <- fread(cache_tsv,sep="\t",header=TRUE,integer64 = "numeric") %>%
          as_tibble %>%
          mutate_at(c("seqnames1","seqnames2"), function(x) factor(x,levels=paste0('chr',c(1:22,"X","Y"))))
        if(!silent) message("Read ",prettyNum(nrow(pixels),big.mark=",")," cached pixels from ",cache_tsv,"...")
      }else{
        pixels  <- cache_tsv
      }
    }
  }

  if(!silent) tic()
  if(!return_table & is.null(cache_tsv)) stop("If a table isn't goint to be returned (return_table), a cache file must be included.")

  if(!is.null(cache_tsv)){
    if(file.exists(cache_tsv) & !overwrite_cache){
      # Note: MAKE SURE fread doesn't do anything with integer64, some of the math starts to chug in R. 'integer64="numeric"' ensures
      # that large integers are turned into doubles instead.
      if(return_table){
        pixels  <- fread(cache_tsv,sep="\t",header=TRUE,integer64 = "numeric") %>%
          as_tibble %>%
          mutate_at(c("seqnames1","seqnames2"), function(x) factor(x,levels=paste0('chr',c(1:22,"X","Y"))))
        if(!silent) message("Read ",prettyNum(nrow(pixels),big.mark=",")," cached pixels from ",cache_tsv,"...")
      }else{
        pixels  <- cache_tsv
      }
    }
  }

  if(is.null(pixels)){
    if(length(file_cool) > 1){
      if(is.null(names(file_cool))){
        names(file_cool)  <- basename(file_cool)
      }

      pixels  <- lapply(names(file_cool),function(fl_nm) {
        pxls  <- read_cooler_hdf5_mc(
          file_cool = file_cool[fl_nm],
          gr_range1=gr_range1,
          gr_range2=gr_range2,
          diag_distance=diag_distance,
          silent=silent,
          max_pixels=max_pixels,
          cache_tsv=NULL) %>%
          mutate(cooler=fl_nm)
      })

      fix_col_names   <- function(pix_tbl,col_names_req=col_nms){
        # Make sure that all retrieved pixel files have all columns (and in the same order).
        # If a table is missing one or more columns, add one with NA values.
        missing_names <- setdiff(col_names_req,colnames(pix_tbl))
        if(length(missing_names) > 0){
          tb_fix <- lapply(missing_names, function(nm){
            tibble(!!as.name(nm):=rep(NA,nrow(pix_tbl)))
          }) %>% do.call(cbind,.)
          pix_tbl_out <- cbind(pix_tbl,tb_fix)
        }else{
          pix_tbl_out <- pix_tbl
        }
        return(pix_tbl_out[,col_names_req])
      }

      col_nms <- lapply(pixels,colnames) %>% unlist %>% unique
      pixels  <- lapply(pixels,fix_col_names,col_names_req = col_nms) %>%
        do.call(rbind,.)
    }
  }
  if(is.null(pixels)){
    # Fetch bin data.
    gr_bins   <- read_cooler_bins_hdf5(file_cooler = file_cool,disable_file_lock = TRUE)

    # Subset relevant bins.
    if(is.null(gr_range1)){
      xbins   <- as.double(gr_bins$bin_id)
      ybins   <- as.double(gr_bins$bin_id)
    }else{
      #If a diagonal distance has been specified, expand gr1 to include 1x that distance on either end, otherwise left and right corners will be omitted in rotated plots.
      if(!is.null(diag_distance)){
        gr1   <- resize(gr_range1,width = width(gr_range1) + 2 * diag_distance,fix = "center")
      }else{
        gr1   <- gr_range1
        #diag_distance <- Inf
      }
      # NOTE: Add one bin to the end of the x-bins list otherwise last row will happen at the start of the last bin.
      xbins<- subsetByOverlaps(gr_bins,gr1,minoverlap = 1L)$bin_id
      xbins<- c(xbins,max(xbins+1))
      if(is.null(gr_range2)){
        ybins <- xbins
      }else{
        ybins <- subsetByOverlaps(gr_bins,gr_range2,minoverlap = 1L)$bin_id
        ybins <- c(ybins,max(ybins + 1))
      }
    }

    # Determine how many pixels are about to be gathered, stop if there are too many.
    pix_num   <- as.double(length(xbins)) * as.double(length(ybins))
    if(pix_num > max_pixels){
      stop("Cooler subset will return up to ",prettyNum(pix_num,big.mark = ",")," pixels...to override this warning and run anyway, increase <max_pixels> or set to 'Inf'.")
    }else{
      if(!silent) message("Fetching up to ",prettyNum(pix_num,big.mark = ",")," pixels...")
    }

    # Load index/bin data.
    hdf5  <- H5Fopen(file_cool,flags="H5F_ACC_RDONLY")
    idx   <- h5read(hdf5,"indexes/bin1_offset") %>%
      tibble(bin_offset = .) %>%
      mutate(bin_offset = as.integer(bin_offset + 1),
             bin_id = as.integer(row_number() - 1))
    tb_bins <- h5read(hdf5,"bins") %>%
      lapply(as.double) %>%
      as_tibble %>%
      mutate(bin_id = row_number())
    H5Fclose(hdf5)

    # Merge bin data into index.
    idx   <- left_join(idx,tb_bins,by="bin_id")

    # Get all bin combinations, and flip those where y < x (bottom triangle) as they are not stored in the .cool matrices.
    tic()
    pixels  <- tibble(binx=xbins,biny=list(xbins)) %>%
      unnest(biny)

    # Get upper-triangle.
    if(is.null(gr_range1) & is.null(gr_range2)){
      # Shortcut for full-genome: don't need to worry about missing lower-triangle pix. Just filter out y<x.
      pixels <- filter(pixels,binx <= biny) %>% mutate(idx = row_number())
    }else{
      # If subset of genome is selected, possible to select lower triangle pix and thus miss them. Must swap x/y coords in this case.
      pixels <- swap_columns(pixels,"binx","biny",function(x,y) x>y) %>% distinct %>% mutate(idx = row_number())
    }
    # Cut pixels into chunks if n_chunks>1.
    if(n_chunks > 1){
      pixels <- pixels %>% mutate(group = cut(idx,breaks = n_chunks,labels = FALSE,include.lowest = TRUE))
    }else{
      pixels <- pixels %>% mutate(group=1)
    }
    pixels <- pixels %>%
      group_by(group) %>%
    # Save chunk results to tempfile() so that table math is kept reasonably small.
      mutate(temp_file = tempfile()) %>%
      ungroup %>%
      group_split(group) %>%
      mclapply(mc.silent=FALSE,mc.cores=n_cores,mc.cleanup = TRUE,FUN = function(tb_in){
        xbins <- tb_in$binx
        ybins <- tb_in$biny

        h5disableFileLocking()
        hdf5  <- H5Fopen(file_cool,flags="H5F_ACC_RDONLY")

        # Determine which rows to retrieve.
        row_range<- filter(idx,bin_id %in% xbins) %>%
          select(bin_offset) %>%
          range

        # Fetch relevant rows in x-dimension then filter y-dimension coordinates.
        pixels  <- h5read(hdf5,name = "pixels",start=row_range[1],count=diff(row_range),bit64conversion="int")

        # Close HDF5.
        H5Fclose(hdf5)

        # Process pixel tibble.
        pixels %>%
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
          mutate(log10_count = log10(count + 1)) %>%
          select(seqnames1,start_adj1,end_adj1,seqnames2,start_adj2,end_adj2,
                 count,log10_count,starts_with("count"),#weight1,weight2,balanced,
                 start1,end1,start2,end2,bin1_id,bin2_id,everything()) %>%
          fwrite(tb_in$temp_file[1])
        return(tb_in$temp_file[1])
      }) %>%
      lapply(function(x) {
        fread(x,integer64="numeric") %>% as_tibble
      }) %>% do.call(rbind,.)
    toc()

    if(!is.null(diag_distance)){
      pixels <- filter(pixels,abs(end2 - end1) <= diag_distance)
    }
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
  if(!return_table){
    rm(pixels)
    pixels <- cache_tsv
  }
  return(pixels)
}
