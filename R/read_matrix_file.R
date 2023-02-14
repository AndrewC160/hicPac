#' @title
#' Read matrix file
#'
#' @description
#' Function reads tabix-indexed matrix files (TSV format) using
#' Rsamtools::scanTabix to filter genomic range coordinates. Separate GRanges
#' can be submitted for asymmetrical matrices (for instance, chr2 on the X axis
#' vs. chr1 on the Y). Matrix files are tabix-indexed based on their first bin
#' coordinate ("bin1"), which is arbitrarily treated here as dimension
#' "X". Files are first filtered using tabix in this dimension, then converted
#' to GenomicRanges and filtered using subsetByOverlaps(). In effect, all
#' coordinates that fall within the X dimension are returned by tabix, then
#' filtered further based on Y coordinates. Theoretically tabix filtering is
#' faster, but the difference seems negligible in most instances...the chief
#' advantage with tabix is not reading an entire matrix file into memory.
#'
#' If both ranges_x and ranges_y are NULL, the entire file is read with
#' fread(). If only one coordinate is provided, the other is set equal by
#' default (making a square matrix). For all/one comparisons, submit a GRanges
#' containing the smaller desired region and another containing all
#' chromosomes.
#'
#' To minimize the size of matrix files, if coordinates (i,j) are included then
#' (j,i) are not as they are functionally equivalent. However, this means that
#' depending on the sorting algorithms used in the generation of matrix files,
#' coordinates may be missing if their equivalents are reflected across X=Y (if
#' you're looking for (i,j) you'll miss it because it's actually (j,i)). To
#' compensate for this, auto_reflect runs the read_matrix_file() function twice
#' to compensate, once normally and once with X and Y ranges switched. This can
#' be disabled, but you should be confident that your selected ranges are all
#' represented. Auto_reflect is not run if an entire matrix file is read (for
#' obvious reasons).
#'
#' Input matrix columns: seqnamesX, startX, endX, seqnamesY, startY, endY, binIDX, binIDY, count.
#'
#' @param matrix_file_in Gzipped and tabix-indexed TSV file of matrix values with nine columns.
#' @param ranges_x List of GRanges to be included on the X axis.
#' @param ranges_y List of GRanges to be included on the Y axis.
#' @param auto_reflect Should function be applied across X=Y as well? Defaults to TRUE.
#'
#' @import dplyr
#' @import data.table
#' @import magrittr
#' @import tibble
#' @import Rsamtools
#' @import GenomicRanges
#' @import clugPac
#'
#' @export

read_matrix_file<- function(matrix_file_in,ranges_x=NULL,ranges_y=NULL,auto_reflect=TRUE){
  select <- dplyr::select
  rename <- dplyr::rename
  col_nms=c("seqnamesX","startX","endX","seqnamesY","startY","endY","binX","binY","count")
  if(is.null(ranges_x) & is.null(ranges_y)){
    #If both dimensions are null, load the entire table (~40 seconds with 66million lines aka 10K resolution)
    tb_mtx  <- fread(matrix_file_in,col.names = col_nms) %>%
      as_tibble
  }else{
    if(is.null(ranges_y)){
      #If Y dimension missing, set it equal to X (square).
      ranges_y <- ranges_x
    }else if(is.null(ranges_x)){
      #If X dimension missing, set it equal to Y (square).
      ranges_x <- ranges_y
    }
    #Assign names to GRanges.
    if(is.null(names(ranges_x))){
      names(ranges_x) <- as.character(ranges_x)
    }
    if(is.null(names(ranges_y))){
      names(ranges_y) <- as.character(ranges_y)
    }
    tb_mtx  <- scanTabix(matrix_file_in,param=ranges_x) %>%
      lapply(function(txt_in){
        if(length(txt_in)==1){
          txt_in <- paste0(txt_in,"\n") #Glitch(?) if one line of text is submitted to fread(text) argument; must have at least one newline.
        }
        fread(text=txt_in,sep="\t",col.names = col_nms) %>%
          as_tibble
      })
    names(tb_mtx)   <- names(ranges_x)
    ranges_y$nms    <- names(ranges_y)

    tb_mtx  <- lapply(names(tb_mtx),function(tb_nm){
      mutate(tb_mtx[[tb_nm]],regionX = tb_nm)
    }) %>%
      do.call(rbind,.) %>%
      mutate_at(c("seqnamesX","seqnamesY"), function(x) factor(x,levels=paste0("chr",c(1:22,"X","Y")))) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                               seqnames.field = "seqnamesY",
                               start.field = "startY",
                               end.field = "endY") %>%
      subsetByOverlaps(ranges_y) %>%
      annotate_gr(gr_subject = ranges_y,cols_query = "regionY",cols_subject = "nms",first_only = TRUE) %>%
      as_tibble %>%
      rename(seqnamesY = seqnames,
             startY = start,
             endY = end) %>%
      select(-width)
    if(auto_reflect){
      tb_mtx <- read_matrix_file(matrix_file_in = matrix_file_in,
                                 ranges_x=ranges_y,
                                 ranges_y=ranges_x,
                                 auto_reflect = FALSE) %>% #Don't auto_reflect recursively.
                  select(colnames(tb_mtx)) %>%
        rbind(tb_mtx,.)
    }
  }
  tb_mtx <- mutate_at(tb_mtx,c("seqnamesX","seqnamesY"),function(x) factor(x,levels=paste0("chr",c(1:22,"X","Y"))))
  return(tb_mtx)
}
