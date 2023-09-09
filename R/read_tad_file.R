#' @title Read TAD file
#' 
#' @description Read a TAD boundary file produced by Cooltools via the hic.cgp
#' pipeline. By default returns a table of boundaries, but can also return a 
#' GRanges object of tads themselves (i.e. the region between two boundaries).
#' Boundary calls are reported as TRUE/FALSE in columns that are named in the 
#' format "x3" to indicate "3 x bin size." For instance, boundary calls measured
#' with in a 100kb resolution matrix with a window size of 300k will be reported
#' in the column "x3," as will a 10kb matrix with a 30kb window size.
#' 
#' @param file_in Filename to read.
#' @param return_tad_gr Should the TAD start/end positions themselves be returned? Defaults to FALSE, in which case each boundary is returned on one row.
#' @param boundary_column Should a specific window size be used to determine boundaries (e.g "x5" for the 5x bin size window, etc)? Defaults to NULL.
#' @param boundary_cols How many boundary calls should be TRUE for a boundary to be returned? Defaults to 1, in which case at least one of the four window sizes must have been defined as a boundary. Can be set to zero, in which case all pixels are returned.
#' 
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import tibble
#'
#' @export

read_tad_file <- function(file_in,return_tad_gr=FALSE,boundary_column=NULL,boundary_calls=1){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  select  <- dplyr::select
  first   <- dplyr::first
  arrange <- dplyr::arrange

  tb_tads <- fread(file_in) %>%
    as_tibble %>%
    filter(!is_bad_bin) %>%
    select(chrom,start,end,starts_with("is_boundary")) %>%
    rename(seqnames=chrom) %>%
    pivot_longer(cols=starts_with("is_boundary"),names_to = "window",names_pattern = "is_boundary_(.+)",values_to="boundary") %>%
    filter(boundary) %>%
    mutate(res = first(end) - first(start),
           window = paste0("x",as.integer(window)/res),
           name = str_extract(basename(file_in),"^[[:alnum:]-]+")) %>%
    pivot_wider(id_cols = c(seqnames,start,end,name,res),names_from=window,values_from=boundary,values_fill = FALSE) %>%
    mutate(total = select(.,starts_with("x")) %>% as.matrix %>% rowSums)
  
  if(!is.null(boundary_column)){
    if(!boundary_column %in% colnames(tb_tads)){
      stop("<boundary_col> must be a column in the tad output file, i.e. 'x5' etc...")
    }
    tb_tads <- filter(tb_tads,!!as.name(boundary_column))
  }else{
    tb_tads <- filter(tb_tads,total >= boundary_calls)
  }
  tb_tads <- tb_tads %>%
    as_tibble %>%
    select(seqnames,start,end,name,res,x3,x5,x10,x25,total)
  
  if(return_tad_gr){
    tb_tads <- tb_tads %>%
      arrange(seqnames,start,end) %>%
      group_by(seqnames) %>%
      mutate(
        pos = (start+end)/2,
        start = pos,
        end = lead(pos),
        start_calls = total,
        end_calls = lead(total)) %>%
      ungroup %>%
      select(seqnames,start,end,name,res,start_calls,end_calls) %>%
      filter(!is.na(end))
  }
  return(tb_tads)
}
