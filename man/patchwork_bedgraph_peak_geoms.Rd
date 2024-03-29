% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchwork_bedgraph_peak_geoms.R
\name{patchwork_bedgraph_peak_geoms}
\alias{patchwork_bedgraph_peak_geoms}
\title{Patchwork bedGraph peak geoms.}
\usage{
patchwork_bedgraph_peak_geoms(
  regions_tb,
  peak_gr,
  epitope_column,
  respect_facets = TRUE,
  peak_alpha = 0.25,
  peak_color_override,
  outline_color = NA,
  outline_linewidth = 0.5,
  outline_linetype = "dashed"
)
}
\arguments{
\item{regions_tb}{Region table, generally the output of patchwork_bins(boundaries_only=TRUE).}

\item{peak_gr}{GRanges object with peaks.}

\item{epitope_column}{Metadata column of peak_gr which differentiates peaks from different epitopes. Typically corresponds to epitope names used in patchwork_bedgraph().}

\item{respect_facets}{Should facetting in patchwork_bedgraph() be respected? Defaults to TRUE, otherwise peaks are drawn over all facets.}

\item{peak_alpha}{Transparency of each peak, defaults to 0.25.}

\item{peak_color_override}{Color to apply to all peaks in the even the scale_fill_manual() values in the patchwork_bedgraph() plot are to be ignored.}

\item{outline_color}{Color of lines to draw around each peak. Defaults to NA, i.e. outlines are not drawn.}

\item{outline_linetype}{Linetype of outlines, if drawn. Defaults to "dashed".}

\item{outline_width}{Linewidth of outlines, if drawn. Defaults to 0.5.}
}
\description{
Returns peak annotations to be applied to patchwork_bedgraph()
plot which respect segment orientations etc. If an "epitope" column is
included, this value will respect the facets and color schemes of the same
column in the patchwork_bedgraph() output plot's aesthetic.
}
