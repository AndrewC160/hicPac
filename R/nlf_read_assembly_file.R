#' @title NeoLoopFinder read assembly file
#' 
#' @description
#' Read in and parse an output assembly file from NeoLoopFinder. Reads file line
#' by line and applies nfl_parse_assembly() to each line, returning a single
#' table with assemblies split up with one segment per line.
#' 
#' @param fl_nm Name of file to read in.
#' 
#' @import dplyr
#' @import tidyr
#' 
#' @export 

nlf_read_assembly_file<- function(fl_nm){
  txt <- scan(fl_nm,what=character(),quiet=TRUE,sep="\n")
  lapply(txt,nlf_parse_assembly) %>%
    do.call(rbind,.)
}