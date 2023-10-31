#' @title Normalize HiC score column
#'
#' @description When plotting multiple HiC interaction matrices together, score
#' can (should?) be normalized based on reads in each matrix. This is best
#' accomplished by taking the sum of score in both samples and dividing each by
#' the lowest sum. Score columns are then divided by their respective ratio.
#' For instance, samples with 1000 and 2000 reads each will have ratios of 1
#' and 2, respectively; score columns are then divided by 1 and 2, ultimately
#' reducing the scores in sample 2 by half.
#'
#' @param tb_pix Pixel tibble to normalize.
#' @param score_column Score column to normalize, defaults to "log10_count_wgs_nrm".
#'
#' @import dplyr
#' @import tibble
#' @import clugPac
#' @import magrittr
#'
#' @export

normalize_score_column  <- function(tb_pix,score_column="log10_count_wgs_nrm"){
  if(!score_column %in% colnames(tb_pix)) stop("Score column '",score_column,"' not found.")
  score_sums  <- tb_pix %>%
    group_by(cooler) %>%
    summarize(score_sum := sum(!!as.name(score_column)),.groups='drop') %>%
    vectify(score_sum,cooler)
  score_ratio <- score_sums / min(score_sums)
  mutate(tb_pix,!!as.name(score_column) := !!as.name(score_column) / score_ratio[cooler]) %>%
    return()
}
