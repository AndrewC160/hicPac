#' @title Cooler to InteractionSet
#'
#' @description
#' Read and convert a Cooler file (or list of files) into an InteractionSet
#' object. Each cooler file will be represented as a single column of each
#' assay.
#'
#' In addition to the default "counts" assay, this function also includes assays
#' for the 'count_wgs_nrm' (WGS-normalized counts) column. Additional assays for
#' filtering (etc.) include:
#'
#' cnA: Copy number for the first bin (anchor).
#' cnB: Copy number for the second anchor.
#' threshold: Expected counts based on same chromosome (distance, logarithmic bins) or different chromosome (median chromosome-chromosome interaction baseline).
#'
#' All assays represent 2-dimensional data (specific to each pixel and sample).
#'
#' Default sample-specific data in colData() includes r_name (sample name),
#' totals (total interactions per sample), and intra_median (inter-chromosome
#' median per sample).
#'
#' Default pixel-specific data in rowData() includes bin1_id and bin2_id, as
#' well as annotation of intra-chromosomal nature (intra=TRUE), raw distance
#' between bin midpoints (distance; NA for inter-chromosomal pixels), and
#' logarithmic distance bins (dist_bin; NA for inter-chromosomal pixels).
#'
#' @param cooler_files Cooler files to be read, one per sample.
#' @param gr_window GenomicRanges object with one or more regions to extract. Defaults to NULL, in which case the entire genome will be read.
#' @param max_pixels Maximum pixels to extract (fed to read_cooler_hdf5); defaults to 1E7.
#' @param col_data_tibble Tibble of column data to merge with the default colData() slot of the IS; defaults to NULL, in which case no data is added.
#' @param col_data_name_field Name field of <col_data_tibble> to merge by; should be the names of each cooler file (which will become column names of assay; see also read_cooler_hdf5() behavior).
#' @param cache_rds RDS filename under which to cache results (reading full genomes can take a very long time). Defaults to NULL, in which case no cache is produced.
#' @param o_write Should existing cache files be over-written? Defaults to FALSE.
#' @param o_write_pixels Should an existing pixel cache file (<cache_rds> but with the ".rds" ending replaced with "_pix.tsv") be over-written? Defaults to FALSE, and no effect if a cache file exists and is not being over-written.
#'
#' @import InteractionSet
#' @import clugPac
#' @import tidyverse
#' @import GenomicRanges
#' @import data.table
#' @import magrittr
#' @import scales
#'
#' @export

cooler_to_InteractionSet  <- function(cooler_files,gr_window = NULL,max_pixels=1e7,col_data_tibble=NULL,col_data_name_field="sample_name",cache_rds=NULL,o_write=FALSE,o_write_pixels=FALSE){
  if(is.null(cache_rds)){
    cache_fl <- ""
  }else{
    cache_fl <- cache_rds
  }
  if(file.exists(cache_fl) & !o_write){
    iset <- readRDS(cache_fl)
  }else{
    tic()
    pix <- read_cooler_hdf5(cooler_files,
                            gr_range1 = gr_window,
                            max_pixels = max_pixels,
                            cache_tsv = gsub(".rds","_pix.tsv",cache_fl),
                            overwrite_cache = o_write_pixels)
    reg <- read_cooler_bins_hdf5(file_cooler = cooler_files[[1]],granges_list = gr_window)
    bin_ids <- reg$bin_id
    regs <- tibble(bin1_id = rep(bin_ids,each=length(bin_ids)),
                   bin2_id = rep(bin_ids,length.out = length(bin_ids)^2)) %>%
      left_join(by="bin1_id",
                as_tibble(reg) %>%
                  select(seqnames,start,end,bin_id) %>%
                  rename_all(paste0,"1") %>%
                  rename(bin1_id=bin_id1)) %>%
      left_join(by="bin2_id",
                as_tibble(reg) %>%
                  select(seqnames,start,end,bin_id) %>%
                  rename_all(paste0,"2") %>%
                  rename(bin2_id=bin_id2)) %>%
      filter(bin2_id >= bin1_id) %>%
      mutate(pix_id = paste0(bin1_id,"_",bin2_id)) %>%
      select(pix_id,bin1_id,bin2_id,everything())

    # Raw counts.
    c1 <- select(pix,bin1_id,bin2_id,cooler,count) %>%
      as_tibble %>%
      pivot_wider(id_cols=c(bin1_id,bin2_id),
                  names_from = cooler,
                  values_from = count,
                  values_fill = 0) %>%
      arrange(bin1_id,bin2_id) %>%
      mutate(pix_id = paste0(bin1_id,"_",bin2_id))

    # Remove bins with no counts across all samples.
    regs <- filter(regs,pix_id %in% c1$pix_id)

    all_regs <- GInteractions(mode="reverse",
                              anchor1=makeGRangesFromDataFrame(regs,seqnames.field = "seqnames1",start.field = "start1",end.field="end1",keep.extra.columns = F),
                              anchor2=makeGRangesFromDataFrame(regs,seqnames.field = "seqnames2",start.field = "start2",end.field="end2",keep.extra.columns = F))

    c1 <- c1 %>%
      select(-bin1_id,-bin2_id,-pix_id) %>%
      as.matrix

    # WGS-normalized counts (from pipeline).
    if("count_wgs_nrm" %in% colnames(pix)){
      c2 <- select(pix,bin1_id,bin2_id,cooler,count_wgs_nrm) %>%
        pivot_wider(id_cols=c(bin1_id,bin2_id),
                    names_from = cooler,
                    values_from = count_wgs_nrm,
                    values_fill = 0) %>%
        arrange(bin1_id,bin2_id) %>%
        select(-bin1_id,-bin2_id) %>%
        as.matrix
    }

    # Pseudo-WGS CNV approximations. For each cooler, sum all counts per bin (column sums, effectively) and divide by the median value (which should be 1).
    # In short, treat HiC reads as single-end reads and count those that fall in each bin (1 dimension).
    bin_cnvs <- select(pix,bin1_id,cooler,count) %>%
      group_by(cooler,bin1_id) %>%
      summarize(count = sum(count),.groups="drop") %>%
      rename(bin_id = bin1_id) %>%
      mutate(cooler = factor(cooler,levels=colnames(c1))) %>%
      group_by(cooler) %>%
      mutate(count = count / median(count)) %>%
      ungroup

    c3a <- regs %>%
      left_join(rename(bin_cnvs,cnA=count),
                by=c("bin1_id"="bin_id"),
                relationship="many-to-many") %>%
      select(cooler,pix_id,bin1_id,cnA) %>%
      arrange(cooler) %>%
      pivot_wider(id_cols = c(pix_id,bin1_id),
                  names_from = cooler,
                  values_from =cnA,
                  values_fill =0) %>%
      select(-pix_id,-bin1_id) %>%
      as.matrix
    if("NA" %in% colnames(c3a)){
      cn  <- which(colnames(c3a) == "NA")
      c3a <- c3a[,-cn]
    }

    c3b <- regs %>%
      left_join(rename(bin_cnvs,cnB=count),
                by=c("bin2_id"="bin_id"),
                relationship="many-to-many") %>%
      select(cooler,pix_id,bin2_id,cnB) %>%
      arrange(cooler) %>%
      pivot_wider(id_cols = c(pix_id,bin2_id),
                  names_from = cooler,
                  values_from =cnB,
                  values_fill =0) %>%
      select(-pix_id,-bin2_id) %>%
      as.matrix
    if("NA" %in% colnames(c3b)){
      cn  <- which(colnames(c3b) == "NA")
      c3b <- c3b[,-cn]
    }

    # CNV-normalized offset values.
    # Define pair-wise bin CNVs by the minimum value of CNVs in anchor A and B as this is the limiting factor of copy number for that pixel.
    # cMin <- pmin(c3a,c3b)

    iset <- InteractionSet(assays = list(counts=c1,wgs_nrm=c2,cnA=c3a,cnB=c3b),interactions = all_regs)

    # Row data.
    ## Intra-chromosomal distance bins.
    max_dist = max(get_seqsizes())
    res <- width(regions(iset)[1])-1
    bin_sizes <- c(0,res)
    while(max(bin_sizes) < max_dist){
      bin_sizes <- c(bin_sizes,bin_sizes[length(bin_sizes)]*2)
    }
    nmsStr <- paste0(prettyBP(bin_sizes),"-")
    nmsStr[length(nmsStr)] <- ">"
    nmsEnd <- prettyBP(bin_sizes)[-1]
    nmsEnd <- c(nmsEnd,nmsEnd[length(nmsEnd)])
    bin_sizes <- setNames(bin_sizes,paste0(nmsStr,nmsEnd))

    rowData(iset) <- cbind(regs,assay(iset,"counts")) %>%
      mutate(intra = seqnames1 == seqnames2) %>%
      rowwise %>%
      mutate(distance = ifelse(intra,(end2+start2)/2 - (end1+start1)/2,NA)) %>%
      ungroup %>%
      mutate(dist_bin = cut(distance,breaks=bin_sizes,labels=FALSE,include.lowest = TRUE,right = FALSE),
             dist_bin = factor(names(bin_sizes[dist_bin]),levels=names(bin_sizes))) %>%
      select(bin1_id,bin2_id,intra,distance,dist_bin)

    # Col data.
    ## Inter-chromosomal cutoffs defined by median values of such bins.
    inter_meds <- colMedians(assay(iset,"counts")[!rowData(iset)$intra,])
    df_col <- colData(iset) %>%
      as.data.frame %>%
      rownames_to_column("r_name") %>%
      mutate(totals = colSums(assay(i = "counts",iset)),
             wgs_nrm_totals = colSums(assay(i="wgs_nrm",iset)),
             intra_median = inter_meds[r_name])

    if(!is.null(col_data_tibble)){
      if(!col_data_name_field %in% colnames(col_data_tibble)) stop(paste0("Name column '",col_data_name_field,"' not present in column data table."))
      df_col <- df_col %>%
        left_join(
          col_data_tibble,
          by=c('r_name'=col_data_name_field)
        ) %>%
        DataFrame
      rownames(df_col) <- df_col$r_name
      colData(iset) <- df_col
    }

    ## Inter/intra-chromosomal cutoff matrix (counts)
    tb_rows <- cbind(anchors(iset,"first") %>% as_tibble %>% select(-strand,-width) %>% rename_all(paste0,"1"),
                     anchors(iset,"second") %>% as_tibble %>% select(-strand,-width) %>% rename_all(paste0,"2"),
                     rowData(iset),
                     assay(iset,"counts")) %>%
      as_tibble

    tb_inter <- tb_rows %>%
      as_tibble %>%
      filter(!intra) %>%
      select(-bin1_id,-bin2_id,-intra,-distance,-dist_bin,-start1,-start2,-end1,-end2) %>%
      pivot_longer(cols = -c(seqnames1,seqnames2),
                   names_to="sample",values_to="count") %>%
      group_by(seqnames1,seqnames2,sample) %>%
      summarize(med = median(count),.groups="drop") %>%
      pivot_wider(id_cols = c(seqnames1,seqnames2),names_from = sample,values_from = med) %>%
      right_join(by=c("seqnames1","seqnames2"),
                 filter(tb_rows,!intra) %>%
                   select(seqnames1,seqnames2,bin1_id,bin2_id)) %>%
      select(bin1_id,bin2_id,everything(),-seqnames1,-seqnames2)

    tb_intra <- tb_rows %>%
      filter(intra) %>%
      select(-seqnames2,-intra,-bin1_id,-bin2_id,-distance,-start1,-end1,-start2,-end2) %>%
      rename(seqnames=seqnames1) %>%
      pivot_longer(cols=-c(seqnames,dist_bin),
                   names_to="sample",values_to="count") %>%
      group_by(sample,seqnames,dist_bin) %>%
      summarize(count = median(count),
                .groups="drop") %>%
      pivot_wider(id_cols = c(seqnames,dist_bin),
                  names_from = sample,values_from=count) %>%
      right_join(by=c("seqnames"="seqnames1","dist_bin"),
                 filter(tb_rows,intra) %>%
                   select(seqnames1,dist_bin,bin1_id,bin2_id)) %>%
      select(bin1_id,bin2_id,everything(),-seqnames,-dist_bin)

    cbck <- rbind(tb_inter,tb_intra) %>%
      arrange(bin1_id,bin2_id) %>%
      select(-bin1_id,-bin2_id) %>%
      as.matrix(dimnames=list(NULL,NULL))

    assay(iset,"threshold") <- cbck

    ## Inter/intra-chromosomal cutoff matrix (WGS-normalized counts).
    if("count_wgs_nrm" %in% colnames(pix)){
      tb_rows <- cbind(anchors(iset,"first") %>% as_tibble %>% select(-strand,-width) %>% rename_all(paste0,"1"),
                       anchors(iset,"second") %>% as_tibble %>% select(-strand,-width) %>% rename_all(paste0,"2"),
                       rowData(iset),
                       assay(iset,"wgs_nrm")) %>%
        as_tibble

      tb_inter <- tb_rows %>%
        as_tibble %>%
        filter(!intra) %>%
        select(-bin1_id,-bin2_id,-intra,-distance,-dist_bin,-start1,-start2,-end1,-end2) %>%
        pivot_longer(cols = -c(seqnames1,seqnames2),
                     names_to="sample",values_to="count") %>%
        group_by(seqnames1,seqnames2,sample) %>%
        summarize(med = median(count),.groups="drop") %>%
        pivot_wider(id_cols = c(seqnames1,seqnames2),names_from = sample,values_from = med) %>%
        right_join(by=c("seqnames1","seqnames2"),
                   filter(tb_rows,!intra) %>%
                     select(seqnames1,seqnames2,bin1_id,bin2_id)) %>%
        select(bin1_id,bin2_id,everything(),-seqnames1,-seqnames2)

      tb_intra <- tb_rows %>%
        filter(intra) %>%
        select(-seqnames2,-intra,-bin1_id,-bin2_id,-distance,-start1,-end1,-start2,-end2) %>%
        rename(seqnames=seqnames1) %>%
        pivot_longer(cols=-c(seqnames,dist_bin),
                     names_to="sample",values_to="count") %>%
        group_by(sample,seqnames,dist_bin) %>%
        summarize(count = median(count),
                  .groups="drop") %>%
        pivot_wider(id_cols = c(seqnames,dist_bin),
                    names_from = sample,values_from=count) %>%
        right_join(by=c("seqnames"="seqnames1","dist_bin"),
                   filter(tb_rows,intra) %>%
                     select(seqnames1,dist_bin,bin1_id,bin2_id)) %>%
        select(bin1_id,bin2_id,everything(),-seqnames,-dist_bin)

      cbck <- rbind(tb_inter,tb_intra) %>%
        arrange(bin1_id,bin2_id) %>%
        select(-bin1_id,-bin2_id) %>%
        as.matrix(dimnames=list(NULL,NULL))

      assay(iset,"threshold_wgsNrm") <- cbck
    }

    if(!is.null(cache_rds)){
      saveRDS(iset,file=cache_fl)
    }
    toc()
  }
  return(iset)
}
