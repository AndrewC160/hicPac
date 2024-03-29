% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchwork_annotations.R
\name{patchwork_annotations}
\alias{patchwork_annotations}
\title{Patchwork annotations.}
\usage{
patchwork_annotations(
  granges_in,
  regions_tb = NULL,
  hic_file = NULL,
  gr_segments = NULL
)
}
\arguments{
\item{granges_in}{GenomicRanges object containing annotations to transform.}

\item{regions_tb}{Table from patchwork_bins(boundaries_only=TRUE) defining the contig in question.}

\item{hic_file}{Cooler file from which to extract bin information if <regions_tb> is NULL; ignored otherwise.}

\item{gr_segments}{GenomicRanges object containing contig segment information to be used if <regions_tb> is NULL; ignored otherwise.}
}
\description{
Given a GenomicRanges object of genome annotations and a
region table from patchwork_bins(boundaries_only=TRUE),subset and clip the
contents of these ranges to identify features which fall within each contig
segment's range. Annotates truncated features with "trunc_start" and
"trunc_end" columns, and while start/end coordinates are provided, X and Y
coordinates for plotting are found in "bin_start" and "bin_end" columns,
which have contig-specific values. These values also respect reversal:
features which fall on segments that have been added in the 3'-5' direction
are reversed ("start" and "end" coordinates run backwards), and strand
information is inverted ("+" strand genes are labeld "-"). If <regions_tb>
is not provided, <hic_file> and <gr_segments> can be substituted to execute
patchwork_bins() within the function.
}
