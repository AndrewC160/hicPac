% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchwork_region_outlines.R
\name{patchwork_region_outlines}
\alias{patchwork_region_outlines}
\title{Patchwork region outlines}
\usage{
patchwork_region_outlines(tb_regs_in)
}
\arguments{
\item{tb_regs_in}{Tibble of patchwork region boundaries output from patchwork_bins(boundaries_only=TRUE).}
}
\description{
Given a boundaries-only regions table from patchwork_bins(), return a table
with coordinates for outlines that can be drawn as polygons around each bin
combination in a rotated HiC interaction matrix. Plot using:

geom_polygon(data=tb_outlines,mapping=aes(x=x,y=y,group=region),inherit.aes=FALSE)

Also, either filter out negative coordinates or use:

scale_y_continuous(oob=scales::oob_keep)
}
