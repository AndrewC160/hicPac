#' @title Patchwork X-axis
#'
#' @description Given a table defining a given patchwork contig (generally the
#' output of patchwork_bins(boundaries_only=TRUE)), define an X-axis and breaks
#' which are more descriptive. Every segment will include at least a start and
#' end break, and junctions between segments will be merged and annotated with
#' both coordinates (chr1:123,456 - chr4:654,321). Segments which are reversed
#' will include a "(-)" prefix indicating as much. Output is a list which
#' includes a title that is a general description of the contig ("6 contigs,
#' 2MB") and a vector of breakpoint labels named by the X-axis value of the
#' contig-specific coordinate system.
#'
#' Example (reversed segment of Chr1 junctioned to a forward segment of Chr2):
#' 1               100                              400
#' (-)chr1:300,000 (-)chr1:100,000-chr2:100,000,000 chr2:100,600,000
#'
#' @param tb_region Table defining the regions of the contig, generally the output of patchwork_bins(boundaries_only=TRUE).
#' @param simplify_labels Should genome coordinates be simplified into basepair units (i.e. 15,000,000 becomes 15MB)? Defaults to FALSE.
#' @param simplify_digits If <simplify_labels> is TRUE, how many digits should be used when rounding? Defaults to 3.
#' @param invert_junction_labels Boolean; should junction labels be reversed? Useful if axis will be used for Y-axis. Defaults to FALSE.
#'
#' @import dplyr
#' @import tidyr
#' @import clugPac
#' @import magrittr
#'
#' @export

patchwork_xaxis   <- function(tb_region,simplify_labels=FALSE,simplify_digits=3,invert_junction_labels=FALSE){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange

  tb  <- tb_region %>%
    filter(!trans_region) %>%
    mutate(width = end1 - start1,
           frac = width / sum(width)) %>%
    select(-seqnames2,-start2,-end2,-bin_alt2_1,-bin_alt2_2,-strand2) %>%
    mutate(frac = width / sum(width))
  x_ttl   <- summarize(tb,count=n(),size=prettyBP(sum(width))) %>% unlist
  x_ttl   <- paste0(x_ttl['count'],
                    ifelse(nrow(tb) > 1," segments, "," segment, "),
                    x_ttl['size'])

  if(simplify_labels){
    lab_func <- function(x) prettyBP(x,digits=simplify_digits)
  }else{
    lab_func <- comma
  }

  x_brks  <- tb %>%
    rowwise %>%
    mutate(x_brks = case_when(frac > 0.7 ~ 5,
                              frac > 0.6 ~ 4,
                              frac > 0.2 ~ 3,
                              TRUE ~ 2)) %>%
    mutate(x_brks = list(as.integer(seq(bin_alt1_1,bin_alt1_2,length.out=x_brks)))) %>%
    unnest(x_brks) %>%
    group_by(region) %>%
    mutate(x_frac = (x_brks - min(x_brks)) / (max(x_brks) - min(x_brks)),
           x_vals = (max(end1) - min(start1))*x_frac + min(start1)) %>%
    group_by(region) %>%
    # 24DEC11: In cases where only one pixel is present per region, x-labels break.
    # In these cases manually establish start/end coordinates/fractions.
    mutate(x_frac = case_when(is.na(x_frac) & row_number() == 1 ~ 0,
                             is.na(x_frac) & row_number() == n() ~ 1,
                             TRUE ~ x_frac),
           x_vals = case_when(is.na(x_vals) & row_number() == 1 ~ start1,
                              is.na(x_vals) & row_number() == n() ~ end1,
                              TRUE ~ x_vals)) %>%
    ##########################################################################
    ungroup %>%
    group_by(region,strand1,seqnames1) %>%
    summarize(x_brks = list(x_brks),
              x_vals = list(x_vals),
              reg_start = min(start1),
              reg_end = max(end1),
              .groups="drop") %>%
    rowwise %>%
    mutate(width = reg_end - reg_start,
           x_vals = ifelse(strand1 == "-",list(rev(x_vals)),list(x_vals))) %>%
    unnest(c(x_brks,x_vals)) %>%
    group_by(region,strand1,seqnames1) %>%
    mutate(x_frac = (x_vals - min(x_vals)) / (max(x_vals) - min(x_vals)),
           x_vals = reg_start + (width * x_frac),
           bound = case_when(row_number() == 1 ~ "left",
                             row_number() == n() ~ "right",
                             TRUE ~ "mid"),
           label = lab_func(x_vals),
           label = case_when(strand1 == "-" & bound != "mid" ~ paste0("(-)",seqnames1,":",label),
                             bound != "mid" ~ paste0(seqnames1,":",label) ,
                             TRUE ~ label)) %>%
    ungroup %>%
    mutate(label = case_when(bound == "right" & lead(bound,default="")=="left" & invert_junction_labels ~ paste(lead(label),label,sep="-\n"),
                             bound == "right" & lead(bound,default="")=="left" ~ paste(label,lead(label),sep="-\n"),
                             TRUE ~ label),
           label = ifelse(bound == "left" & lag(bound,default="")=="right","",label)) %>%
    filter(label != "") %>%
    vectify(x_brks,label)

  x_lims <- range(x_brks) + c(0,1)
  return(list(title = x_ttl,breaks=x_brks,limits=x_lims))
}

