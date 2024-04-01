#' @title InteractionSet to pixels
#' 
#' @description Given an InteractionSet object, convert it into a pixel table 
#' compaitible with other hicPac features. All assays will be slotted into one
#' column of a long-format pixel table, and all assay columns will be assigned 
#' cooler names in the "cooler" column.
#' 
#' @param iset_in InteractionSet object to convert.
#' 
#' @import diffHic
#' @import InteractionSet
#' @import tidyr
#' @import magrittr
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import tibble
#' @import clugPac
#' 
#' @export

InteractionSet_to_pixels <- function(iset_in){
  seq_adj <- get_seqsizes_adj()
  tb_out  <- cbind(
    anchors(iset_in,"first") %>% as_tibble() %>% select(-strand,-width) %>% rename_all(paste0,"1"),
    anchors(iset_in,"second") %>% as_tibble() %>% select(-strand,-width) %>% rename_all(paste0,"2"),
    rowData(iset_in)
  ) %>%
    as_tibble %>%
    mutate(start_adj1 = start1 + seq_adj[as.character(seqnames1)],
           end_adj1 = end1 + seq_adj[as.character(seqnames1)],
           start_adj2 = start2 + seq_adj[as.character(seqnames2)],
           end_adj2 = end2 + seq_adj[as.character(seqnames2)])
  
  # Read all assays into separate long-format tibbles.
  tbs_assays  <- lapply(names(assays(iset_in)), function(as_nm){
    mtx <- cbind(rowData(iset_in)[,c("bin1_id","bin2_id")],
                 assay(iset_in,as_nm)) %>%
      as_tibble
    colnames(mtx) <- c("bin1_id","bin2_id",colData(iset_in)$r_name)
    pivot_longer(mtx,
                 cols = -c(bin1_id,bin2_id),
                 names_to = "cooler",
                 values_to = as_nm)
  })
  # Bind tibbles using cbind (don't join since it'll take a lot longer and assay
  # tables should all be identical, dimension-wise).
  tb_assays <- tbs_assays[[1]]
  for(asy in tbs_assays[-1]){
    tb_assays <- cbind(tb_assays,select(asy,-bin1_id,-bin2_id,-cooler))
  }
  
  # Select basic columns to maintain order consistent with default pixel tables,
  # but let counts and all non-standard assays fall where they may.
  tb_out <- left_join(tb_out,tb_assays,by=c("bin1_id","bin2_id")) %>%
    select(seqnames1,start_adj1,end_adj1,seqnames2,start_adj2,end_adj2,start1,end1,start2,end2,bin1_id,bin2_id,everything())
  return(tb_out)
}