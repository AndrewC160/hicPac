#' @title Rotate patchwork annotations 45 degrees.
#'
#' @description Given a table of annotations specific to a patchwork (generally
#' the output of patchwork_annotations()), calculate "rotated" version which
#' include X1 and X2 coordinates and a Y coordinate denoting their intersection
#' in a rotated HiC matrix. Appends columns x1 (bin of first annotation), x2
#' (mid point between first and second annotation), x3 (bin of second
#' annotation, and y2 (point on Y-axis where junction should reside; y1 and y3
#' in this context would both be zero). Annotations can have a name column
#' (specified via <name_field>), but otherwise they will be assigned names in
#' the format of "annot_1", "annot_2", etc.
#'
#' NOTE: IDs may not be unique in cases where the same genome segment is
#' included multiple times; these are allowed with this function, but output
#' tables may be confusing when designing plots if you forget this detail.
#'
#' @param features_tb Table of genomic features which has at minimum two seqnames and two start fields.
#' @param regions_tb Regions table (output of hicPac::patchwork_bins(boundaries_only=TRUE)).
#' @param gr_list Contig segment list from which to retrieve segment coordinates if <regions_tb> is not provided.
#' @param hic_file Cooler file from which to retrieve bin data if <regions_tb> is not provided.
#' @param seqnames1_field Name of column for chromosome in first annotation, defaults to "seqnames1".
#' @param start1_field Name of column for BP coordinate in first annotation, defaults to "start1".
#' @param seqnames2_field Name of column for chromosome in second annotation, defaults to "seqnames2".
#' @param start2_field Name of column for BP coordinate in second annotation, defaults to "start2".
#' @param name_field Optional; name of column containing name/id data that will be carried over into the output. If not provided, annotations will be labeled "annot1" etc.
#' @param id_fields Optional; names of columns containing ID data that is identical for both variants and should be included in the output table (cell line, SV type, etc.).
#'
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import magrittr
#' @import clugPac
#'
#' @export

patchwork_rotate_annotations  <- function(features_tb,regions_tb=NULL,gr_list=NULL,hic_file=NULL,seqnames1_field="seqnames1",start1_field="start1",seqnames2_field="seqnames2",start2_field="start2",name_field,id_fields){
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  select  <- dplyr::select
  arrange <- dplyr::arrange

  if(is.null(regions_tb)){
    if(is.null(hic_file) | is.null(gr_list)) stop("Either a regions table from patchwork_bins() is required or both the hic_file and gr_segments arguments must be provided.")
    tb_regs <- patchwork_bins(gr_list = gr_list,hic_file = hic_file)
  }else{
    tb_regs <- regions_tb
  }

  region_grange <- tb_regs %>%
    filter(!trans_region) %>%
    select(seqnames1,start1,end1,strand1,region,bin_alt1_1,bin_alt1_2) %>%
    makeGRangesFromDataFrame(seqnames.field="seqnames1",start.field="start1",end.field="end1",strand.field = "strand1",keep.extra.columns = TRUE)
  if(missing(name_field)){
    features_tb <- mutate(features_tb,sv_name=paste0("annot_",row_number()))
  }else{
    features_tb <- mutate(features_tb,sv_name = !!as.name(name_field))
  }
  if(missing(id_fields)){
    id_fields   <- NULL
  }
  sv1     <- makeGRangesFromDataFrame(features_tb,keep.extra.columns = TRUE,
                                      seqnames.field=seqnames1_field,
                                      start.field=start1_field,
                                      end.field=start1_field)[,c("sv_name",id_fields)]
  sv2     <- makeGRangesFromDataFrame(features_tb,keep.extra.columns = TRUE,
                                      seqnames.field=seqnames2_field,
                                      start.field=start2_field,
                                      end.field=start2_field)[,"sv_name"]
  olaps1  <- findOverlaps(query = sv1,subject = region_grange,ignore.strand=TRUE)
  olaps2  <- findOverlaps(query = sv2,subject = region_grange,ignore.strand=TRUE)

  sv1 <- cbind(
    sv1[queryHits(olaps1)] %>% as_tibble %>% select(-width) %>% rename_at(c("seqnames","start","end","strand"),function(x) paste0("sv_",x,"1")),
    as_tibble(region_grange[subjectHits(olaps1)]) %>%
      select(seqnames,start,end,strand,region,bin_alt1_1,bin_alt1_2)
    )  %>%
    as_tibble %>%
    mutate(strand = ifelse(as.character(strand) == "-","-","+") %>% factor,
           start_frac = (sv_start1 - start)/(end - start),
           start_frac = ifelse(strand == "-",1-start_frac,start_frac),
           bin_alt1 = bin_alt1_1 + (bin_alt1_2 - bin_alt1_1)*start_frac) %>%
    select_at(c("sv_seqnames1","sv_start1","sv_name","bin_alt1",id_fields))

  sv2 <- cbind(
    sv2[queryHits(olaps2)] %>% as_tibble %>% select(-width) %>% rename_at(c("seqnames","start","end","strand"),function(x) paste0("sv_",x,"2")),
    as_tibble(region_grange[subjectHits(olaps2)]) %>%
      select(seqnames,start,end,strand,region,bin_alt1_1,bin_alt1_2)
    )  %>%
    as_tibble %>%
    mutate(strand = ifelse(as.character(strand) == "-","-","+") %>% factor,
           start_frac = (sv_start2 - start)/(end - start),
           start_frac = ifelse(strand == "-",1-start_frac,start_frac),
           bin_alt2 = bin_alt1_1 + (bin_alt1_2 - bin_alt1_1)*start_frac) %>%
    select(sv_seqnames2,sv_start2,sv_name,bin_alt2)

  #Multiple="all" because SVs may fall on segments that are used repeatedly; in this case all are important.
  inner_join(sv1,sv2,by="sv_name",multiple="all") %>%
    select(sv_name,everything()) %>%
    mutate(x1 = bin_alt1,
           x2 = (bin_alt1 + bin_alt2)/2,
           x3 = bin_alt2,
           y2 = abs(bin_alt2 - bin_alt1)/2) %>%
    return()
}
