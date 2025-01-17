% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchwork_xaxis.R
\name{patchwork_xaxis}
\alias{patchwork_xaxis}
\title{Patchwork X-axis}
\usage{
patchwork_xaxis(
  tb_region,
  simplify_labels = FALSE,
  simplify_digits = 3,
  invert_junction_labels = FALSE
)
}
\arguments{
\item{tb_region}{Table defining the regions of the contig, generally the output of patchwork_bins(boundaries_only=TRUE).}

\item{simplify_labels}{Should genome coordinates be simplified into basepair units (i.e. 15,000,000 becomes 15MB)? Defaults to FALSE.}

\item{simplify_digits}{If <simplify_labels> is TRUE, how many digits should be used when rounding? Defaults to 3.}

\item{invert_junction_labels}{Boolean; should junction labels be reversed? Useful if axis will be used for Y-axis. Defaults to FALSE.}
}
\description{
Given a table defining a given patchwork contig (generally the
output of patchwork_bins(boundaries_only=TRUE)), define an X-axis and breaks
which are more descriptive. Every segment will include at least a start and
end break, and junctions between segments will be merged and annotated with
both coordinates (chr1:123,456 - chr4:654,321). Segments which are reversed
will include a "(-)" prefix indicating as much. Output is a list which
includes a title that is a general description of the contig ("6 contigs,
2MB") and a vector of breakpoint labels named by the X-axis value of the
contig-specific coordinate system.

Example (reversed segment of Chr1 junctioned to a forward segment of Chr2):
1               100                              400
(-)chr1:300,000 (-)chr1:100,000-chr2:100,000,000 chr2:100,600,000
}
