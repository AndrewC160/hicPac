#'@title
#'NeoLoopFinder parse assembly
#'
#'@description
#'Given a line of text from a NeoLoopFinder assembly file, parse it into a table
#'which is organized as one line per segment with the assembly ID listed in an 
#'"assembly" column. Strand orientation is represented in the "strand" column, 
#'and each segment is labeled as "seg_1" etc. Columns "type_start" and "type_end"
#'indicate which type of SV produced each junction in the assembly (start on the
#'left and end on the right), and in the case of the left- or right-most ends of
#'the assembly these will be labeled as "left_bound" or "right_bound",
#'respectively.
#'
#'@param txt_in Text line from assembly text file.
#'
#'@import stringr
#'@import dplyr
#'@import magrittr
#'@import tidyr
#'@import tibble
#'
#'@export

nlf_parse_assembly    <- function(txt_in){
  rename<- dplyr::rename
  mutate<- dplyr::mutate
  filter<- dplyr::filter
  select<- dplyr::select
  
  txt   <- str_match(txt_in,"^([:alnum:]+)\t(.+$)")
  as_nm <- txt[2] 
  txt   <- txt[3] %>% str_split("\t",simplify = TRUE)
  rgt_bnd <- txt[length(txt)] %>% str_match("([:alnum:]+),([:digit:]+)$") %>% .[,-1]
  lft_bnd <- txt[length(txt)-1] %>% str_match("([:alnum:]+),([:digit:]+)$") %>% .[,-1]
  txt   <- txt[-c(length(txt),length(txt)-1)]
  types <- str_extract(txt,"^[:alpha:]+") %>% na.omit %>% as.character %>% paste(c(1:length(.)),sep="_")
  tb <- paste(txt,collapse=",") %>%
    str_extract_all("[:alnum:]+,[:digit:]+,[\\+\\-]") %>%
    unlist %>%
    paste(collapse="\n") %>%
    read.table(text = .,sep=",",col.names=c("seqnames","pos","orient")) %>%
    as_tibble %>%
    mutate(joint = rep(c("end","start"),length.out = n()),
           idx = row_number(),
           type = rep(types,each=2))
  
  rbind(
    tibble(seqnames=lft_bnd[1],
           pos=lft_bnd[2],
           orient=tb$orient[1],
           joint="start",
           idx=0,
           type="left_bound"),
    tb,
    tibble(seqnames=rgt_bnd[1],
           pos=rgt_bnd[2],
           orient=tb$orient[nrow(tb)],
           joint="end",
           idx=nrow(tb)+1,
           type="right_bound")) %>%
    mutate(pos = as.double(pos),
           seg = paste0("seg_",rep(c(1:(n()/2)),each=2))) %>%
    pivot_wider(id_cols = seg,
                names_from = joint,
                values_from = c(seqnames,pos,type)) %>%
    rename(seqnames=seqnames_start,
           start=pos_start,
           end=pos_end) %>%
    select(seqnames,start,end,seg,type_start,type_end) %>%
    mutate(seqnames=factor(paste0("chr",seqnames),levels=paste0("chr",c(1:22,"X","Y"))),
           strand = ifelse(end > start,"+","-"),
           assembly = as_nm,
           start_adj = ifelse(start>end,end,start),
           end_adj = ifelse(start>end,start,end),
           start=start_adj,
           end=end_adj) %>%
    select(-start_adj,-end_adj) %>%
    # group_split(seg) %>%
    # lapply(makeGRangesFromDataFrame,keep.extra.columns = TRUE) %>%
    # setNames(sapply(.,function(x) x$seg[1])) %>%
    return()
}