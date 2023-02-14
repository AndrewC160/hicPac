#' @title  HiC translocation matrix
#'
#' @description Retrive and combine three pixel tables describing a point
#' translocation: pixels from the origin sequence (side A), the destination
#' sequence (side B), and the pixels representing their junction (side AB).
#' Pixels are cached similarly to read_cooler_hdf5(), but include extra columns
#' "side," "x_bin,", and "y_bin." Side indicates which section of the genome
#' each pixel originates from, and bin columns represent arbitrary bin IDs
#' relative to the translocation center. If "rotate" is TRUE, this combined
#' table is rotated relative to the left-most coordinate based on "x_bin" and
#' "y_bin" coordinates.
#'
#' If the argument 'window_size' is provided, GRange inputs are used only to
#' specify the chromosome and translocation positions, then windows of the same
#' size are generated as appropriate. If 'window_fix' is set as "center" (the
#' default value), the center point of each GRange is expanded to the specified
#' window size (resize(fix="center")). If set to "junction", the start/end
#' locations of the genomic coordinates are assumed to represent the "exact"
#' position of the translocation, and windows are expanded in accordance with
#' that and with orientation. For instance, if chromosome A and B have a
#' translocation at 1,000,000 and 2,000,000, respectively, and both are
#' + oriented, the translocation point from A and B are the End and Start
#' positions, respectively. If A is reversed and B is forward, these positions
#' are the Start and Start positions, respectively (Etc).
#'
#' @param gr_regionA Genomic region to plot on side A.
#' @param gr_regionB Genomic region to plot on side B.
#' @param reverseA If TRUE, A bins will be flipped to place translocation at the right-most point.
#' @param reverseB If TRUE, B bins will be flipped to place translocation at the left-most point.
#' @param window_size Defaults to NULL; if provided, GRanges are configured to the specified window size using the start/end position of input GR regions as appropriate with direction taken into account.
#' @param window_fix Defaults to 'center', in which case input ranges are resized from the center point of the range. If set to "junction," points are taken as the precise point of the translocation, and windows are sized relative to this position and the translocation's orientation.
#' @param cooler_file Cooler file from which to retrieve data.
#' @param rotate Should the output table be rotated 45 degrees? Similar to rot_pix_45(), but using arbitrary x_bin and y_bin columns.
#' @param cache_tsv Cache TSV file for output table.
#' @param o_write Should an existing cache TSV file be overwritten? Defaults to FALSE, in which case the existing table is read.
#'
#' @import tibble
#' @import dplyr
#' @import GenomicRanges
#' @import clugPac
#' @import data.table
#' @import magrittr
#'
#' @export

translocation_matrix<- function(gr_regionA,gr_regionB,reverseA=FALSE,reverseB=FALSE,window_size=NULL,window_fix="center",cooler_file,rotate=TRUE,cache_tsv,o_write=FALSE){
  mutate  <- dplyr::mutate
  rename  <- dplyr::rename
  filter  <- dplyr::filter
  arrange <- dplyr::arrange
  unnest  <- tidyr::unnest

  if(missing(cache_tsv)) cache_tsv <- ""
  if(file.exists(cache_tsv) & !o_write){
    tb_pix    <- fread(cache_tsv,sep="\t",header=TRUE) %>%
      as_tibble() %>%
      mutate_at(c("seqnames1","seqnames2"),function(x) factor(x,levels=paste0("chr",c(1:22,"X","Y"))))
  }else{
    #If window size is provided, ignore GRange sizes and figure out from start/end position.
    if(!is.null(window_size)){
      if(window_fix == "center"){
        gr_regionA  <- resize(gr_regionA,width=window_size,fix="center")
        gr_regionB  <- resize(gr_regionB,width=window_size,fix="center")
      }else if(window_fix == "junction"){
        if(reverseA){
          gr_regionA  <- resize(gr_regionA,width=window_size,fix="start")
        }else{
          gr_regionA  <- resize(gr_regionA,width=window_size,fix="end")
        }
        if(reverseB){
          gr_regionB  <- resize(gr_regionB,width=window_size,fix="start")
        }else{
          gr_regionB  <- resize(gr_regionB,width=window_size,fix="end")
        }
      }
    }
    tb_pix_A  <- read_cooler_hdf5(cooler_file,gr_range1 = gr_regionA) %>% mutate(side = "A")
    tb_pix_B  <- read_cooler_hdf5(cooler_file,gr_range1 = gr_regionB) %>% mutate(side = "B")
    tb_pix_AB <- read_cooler_hdf5(cooler_file,gr_range1 = gr_regionA,gr_range2 = gr_regionB) %>% mutate(side = "AB")
    tb_pix    <- rbind(tb_pix_A,tb_pix_B,tb_pix_AB)

    tb_bins <- read_cooler_bins_hdf5(cooler_file)
    bins_A  <- tb_bins %>%
      subsetByOverlaps(gr_regionA) %>%
      as_tibble %>%
      select(seqnames,start,end,start_adj,end_adj,bin_id) %>%
      mutate(side = "A")
    bins_B  <- tb_bins %>%
      subsetByOverlaps(gr_regionB) %>%
      as_tibble %>%
      select(seqnames,start,end,start_adj,end_adj,bin_id) %>%
      mutate(side = "B")
    #Flip X coordinates if needed.
    if(reverseA){
      bins_A  <- mutate(bins_A,bin=-1*(bin_id - min(bin_id) - 0.5)) %>%
        arrange(bin)
    }else{
      bins_A  <- mutate(bins_A,bin=bin_id - max(bin_id) - 0.5)
    }
    if(reverseB){
      bins_B  <- mutate(bins_B,bin=-1*(bin_id - max(bin_id) - 0.5)) %>%
        arrange(bin)
    }else{
      bins_B  <- mutate(bins_B,bin=bin_id - min(bin_id) + 0.5)
    }

    tb_bins <- rbind(bins_A,bins_B)
    x_rng   <- c(min(tb_bins$bin),max(tb_bins$bin))
    y_rng   <- c(0,max(tb_bins$bin))

    tb_pix <- tb_pix %>%
      inner_join(
        select(tb_bins,bin_id,bin),
        by=c("bin1_id"="bin_id")
      ) %>%
      rename(x_bin = bin) %>%
      inner_join(
        select(tb_bins,bin_id,bin),
        by=c("bin2_id"="bin_id")) %>%
      rename(y_bin = bin) %>%
      #If coordinates have flipped over the diagonal, switch X and Y.
      mutate(x_bin2=ifelse(x_bin > y_bin,y_bin,x_bin),
             y_bin2=ifelse(x_bin > y_bin,x_bin,y_bin),
             x_bin = x_bin2,
             y_bin = y_bin2) %>%
      select(-x_bin2,-y_bin2)
    if(!cache_tsv == ""){
      fwrite(tb_pix,file = cache_tsv,quote = FALSE,row.names=FALSE,col.names = TRUE,sep="\t")
    }
  }
  if(rotate){
    degree  <- pi * -45 / 180
    res     <- 0.5
    origin  <- min(tb_pix$x_bin)

    tb_pix  <- tb_pix %>%
      mutate(x_bin = x_bin - origin,
             y_bin = y_bin - origin,
             x_rot = (x_bin * cos(degree) - y_bin * sin(degree)),
             y_rot = (y_bin * cos(degree) + x_bin * sin(degree)),
             x_rot = x_rot / sqrt(2) + origin,
             y_rot = y_rot / sqrt(2) + origin,
             pix_id = paste(x_bin,y_bin,sep="_")) %>%
      #Pixel coordinate order: Top,Right,Bottom,Left
      mutate(x1=x_rot,
             x2=x_rot + res,
             x3=x_rot,
             x4=x_rot - res,
             y1=y_rot + res,
             y2=y_rot,
             y3=y_rot - res,
             y4=y_rot) %>%
      pivot_longer(cols = c("x1","x2","x3","x4","y1","y2","y3","y4"),
                   names_to = c(".value","num"),
                   names_pattern = "([xy])([1234])") %>%
      select(-num,-x_rot,-y_rot) %>%
      rename(x_coords = x,y_coords=y) %>%
      mutate(y_coords = y_coords - min(y_coords))
  }
  return(tb_pix)
}
