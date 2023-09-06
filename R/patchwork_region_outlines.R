#' @title Patchwork region outlines
#'
#' @description
#' Given a boundaries-only regions table from patchwork_bins(), return a table
#' with coordinates for outlines that can be drawn as polygons around each bin
#' combination in a rotated HiC interaction matrix. Plot using:
#'
#' geom_polygon(data=tb_outlines,mapping=aes(x=x,y=y,group=region),inherit.aes=FALSE)
#'
#' Also, either filter out negative coordinates or use:
#'
#' scale_y_continuous(oob=scales::oob_keep)
#'
#' @param tb_regs_in Tibble of patchwork region boundaries output from patchwork_bins(boundaries_only=TRUE).
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#'
#' @export

patchwork_region_outlines <- function(tb_regs_in){
  mutate  <- dplyr::mutate

  tb_regs_in %>%
    mutate(
      widA = bin_alt1_2 - bin_alt1_1,
      widB = bin_alt2_1 - bin_alt1_2,
      widC = bin_alt2_2 - bin_alt2_1,
      xA = (bin_alt1_1 + (widA + widB + widC)/2)-0.5,
      xB = (bin_alt1_1 + widA + (widB + widC)/2),
      xC = (bin_alt1_2 + (widB)/2)-0.5,
      xD = (bin_alt1_1 + (widA + widB)/2)-1,
      yA = ((widA + widB + widC)/2)+0.5,
      yB = ((widB + widC)/2),
      yC = (widB/2)-0.5,
      yD = ((widA+widB)/2)) %>%
    pivot_longer(cols=c(xA,xB,xC,xD,yA,yB,yC,yD),
                 names_pattern=c("(.)(.)"),
                 names_to=c("dim","coord"),
                 values_to = "bin_alt") %>%
    pivot_wider(id_cols = c(region,trans_region,coord),
                values_from = bin_alt,
                names_from = dim) %>%
    group_by(region) %>%
    mutate(midX=mean(x),
           midY=mean(y)) %>%
    ungroup %>%
    return()
}
