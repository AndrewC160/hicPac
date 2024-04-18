#' @title Pack Cooler HDF5
#'
#' @description
#' Given a table with two-dimensional genome coordinates (specifically columns
#' for seqnames1, start1, end1, seqnames2, start2, end2), pack this into a
#' Cooler file. This entails generating bin IDs and an index, but note that
#' other than chromosome sizes these values are generated based on the table
#' input: as such it is not guaranteed that bin IDs will match between coolers.
#'
#' Coordinate columns will be incorporated into bins/indices etc., and all other
#' columns are treated as columns for the pixel group. If the specified cooler
#' file already exists that filename is returned with a warning unless
#' overwrite_hdf5=TRUE, in which case that file will be deleted and a new one
#' generated in its place.
#'
#' @param pixel_tibble Table of data to convert.
#' @param cooler_file Cooler filename to create.
#' @param seq_sizes Named vector of chromosome sizes to use. If missing, hg38 sizes returned by clugPac::get_seqsizes() will be used.
#' @param overwrite_hdf5 Should existing cooler files be destroyed? Defaults to FALSE.
#' @param chunk_size Size of chunks that pixel values should be split into. Defaults to 1E7.
#'
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import rhdf5
#' @import rhdf5filters
#' @import Rhdf5lib
#' @import clugPac
#'
#' @export

pack_cooler_hdf5<- function(pixel_tibble,cooler_file,seq_sizes,overwrite_hdf5=FALSE,chunk_size=1E7){
  rename <- dplyr::rename
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange
  select <- dplyr::select

  if(file.exists(cooler_file)){
    if(overwrite_hdf5){
      file.remove(cooler_file)
    }else{
      message("Returning existing Cooler file '",cooler_file,"'.")
      cl_file <- cooler_file
    }
  }
  if(!file.exists(cooler_file)){
    message("Writing ",cooler_file,"...")
    tb  <- pixel_tibble
    if(missing(seq_sizes)) seq_sizes <- clugPac::get_seqsizes()

    # Check that minimum column names are present.
    c_nms   <- setdiff(
      c("seqnames1","start1","end1","seqnames2","start2","end2"),
      colnames(tb))
    if(length(c_nms) > 0) stop("Table is missing columns for ",oxford_collapse(c_nms),".")

    # Bins.
    tb_bins <- rbind(
      select(tb,seqnames1,start1,end1) %>% rename_all(function(x) gsub("1$","",x)),
      select(tb,seqnames2,start2,end2) %>% rename_all(function(x) gsub("2$","",x))
    ) %>% distinct %>%
      arrange(seqnames,start) %>%
      mutate(bin = row_number())

    # Pixels.
    tb <- tb %>%
      left_join(
        rename_all(tb_bins,paste0,"1")
      ) %>%
      left_join(
        rename_all(tb_bins,paste0,"2")
      ) %>% suppressMessages %>%
      rename(bin1_id=bin1,bin2_id=bin2) %>%
      # Keep seqnames1 until indices have been built.
      select(-start1,-end1,-seqnames2,-start2,-end2)

    # Index offsets (bin1 and chrom).
    b_offs <- tb %>%
      mutate(row_n = row_number()) %>%
      group_by(bin1_id) %>%
      summarize(offset = min(row_n) - 1,.groups="drop") %>%
      vectify(offset,bin1_id)

    c_offs <- tb %>%
      mutate(row_n = row_number()) %>%
      group_by(seqnames1) %>%
      summarize(offset = min(row_n) - 1,.groups="drop") %>%
      vectify(offset,seqnames1)

    # Finalize pixels.
    tb <- select(tb,-seqnames1)

    h5_name <- cooler_file
    hdf5    <- H5Fcreate(h5_name)

    # Write bins.
    h5_bins <- H5Gcreate(hdf5,"bins")
    h5_bins$chrom <- as.character(tb_bins$seqnames)
    h5_bins$start <- as.double(tb_bins$start)
    h5_bins$end   <- as.double(tb_bins$end)

    # Write chromosomes.
    h5_chrs <- H5Gcreate(hdf5,"chroms")
    h5_chrs$length<- setNames(seq_sizes,NULL)
    h5_chrs$name  <- names(seq_sizes)

    # Write index.
    h5_idx  <- H5Gcreate(hdf5,"indexes")
    h5_idx$bin1_offset <- setNames(b_offs,NULL)
    h5_idx$chrom_offset<- setNames(c_offs,NULL)

    # Write pixels.
    #h5createGroup(hdf5,"pixels")
    h5_pix  <- H5Gcreate(hdf5,"pixels")
    #if(nrow(tb) < chunk_size) chunk_size <- NULL
    for(nm in rev(colnames(tb))){
      v_vals <- tb[[nm]]
      if(!is.numeric(v_vals)){
        v_vals  <- as.character(v_vals)
      }
      h5createDataset(h5_pix,dataset=nm,dims=length(v_vals),chunk=chunk_size)
      h5writeDataset(v_vals,h5_pix,name=nm)
    }
    H5Fclose(hdf5)
  }
  return(cooler_file)
}

