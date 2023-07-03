#' @title NeoLoopFinder read loop file
#' 
#' @description
#' Given a Neoloop text file output by NeoLoopFinder, parse the file into a
#' tibble with one row per loop-assembly combination (the same loop can be 
#' found in multiple lines if it is found in multiple assemblies).
#' 
#' @param file_in NeoLoopFinder Neoloop output file
#' 
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import tibble
#' @import data.table
#' 
#' @export

nlf_read_loop_file    <- function(file_in=fl_neoloops){
  rename<- dplyr::rename
  mutate<- dplyr::mutate
  filter<- dplyr::filter
  select<- dplyr::select
  
  fread(file_in,col.names = c("seqnames1","start1","end1","seqnames2","start2","end2","assemblies")) %>%
    as_tibble %>%
    mutate_at(c("seqnames1","seqnames2"), function(x) factor(x,levels=paste0("chr",c(1:22,"X","Y")))) %>%
    arrange(seqnames1,start1,seqnames2,start2) %>%
    mutate(loop_id = paste0("loop_",row_number())) %>%
    rowwise %>%
    mutate(assemblies = str_split(assemblies,pattern=","),
           assemb_col = list(rep(c("assembly","dist","neo"),length(assemblies)/3))) %>%
    ungroup %>%
    unnest(c(assemblies,assemb_col)) %>%
    pivot_wider(names_from = assemb_col,
                values_from = assemblies,
                values_fn=function(x) paste0(x,collapse=",")) %>%
    rowwise %>%
    mutate_at(c("assembly","dist","neo"), function(x) str_split(x,pattern=",",simplify=FALSE)) %>%
    unnest(c(assembly,dist,neo)) %>%
    mutate(neo = as.logical(as.integer(neo)),
           res = end1 - start1,
           pos1 = (start1 + end1)/2,
           pos2 = (start2 + end2)/2) %>%
    select(seqnames1,pos1,seqnames2,pos2,loop_id,assembly,dist,neo,res) %>%
    return()
}