% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translocation_matrix.R,
%   R/translocation_matrix_v1.R
\name{translocation_matrix}
\alias{translocation_matrix}
\title{HiC translocation matrix}
\usage{
translocation_matrix(
  gr_regionA,
  gr_regionB,
  reverseA = FALSE,
  reverseB = FALSE,
  window_size = NULL,
  window_fix = "center",
  cooler_file,
  rotate = TRUE,
  cache_tsv,
  o_write = FALSE
)

translocation_matrix(
  gr_regionA,
  gr_regionB,
  reverseA = FALSE,
  reverseB = FALSE,
  window_size = NULL,
  window_fix = "center",
  cooler_file,
  rotate = TRUE,
  cache_tsv,
  o_write = FALSE
)
}
\arguments{
\item{gr_regionA}{Genomic region to plot on side A.}

\item{gr_regionB}{Genomic region to plot on side B.}

\item{reverseA}{If TRUE, A bins will be flipped to place translocation at the right-most point.}

\item{reverseB}{If TRUE, B bins will be flipped to place translocation at the left-most point.}

\item{window_size}{Defaults to NULL; if provided, GRanges are configured to the specified window size using the start/end position of input GR regions as appropriate with direction taken into account.}

\item{window_fix}{Defaults to 'center', in which case input ranges are resized from the center point of the range. If set to "junction," points are taken as the precise point of the translocation, and windows are sized relative to this position and the translocation's orientation.}

\item{cooler_file}{Cooler file from which to retrieve data.}

\item{rotate}{Should the output table be rotated 45 degrees? Similar to rot_pix_45(), but using arbitrary x_bin and y_bin columns.}

\item{cache_tsv}{Cache TSV file for output table.}

\item{o_write}{Should an existing cache TSV file be overwritten? Defaults to FALSE, in which case the existing table is read.}
}
\description{
Retrive and combine three pixel tables describing a point
translocation: pixels from the origin sequence (side A), the destination
sequence (side B), and the pixels representing their junction (side AB).
Pixels are cached similarly to read_cooler_hdf5(), but include extra columns
"side," "x_bin,", and "y_bin." Side indicates which section of the genome
each pixel originates from, and bin columns represent arbitrary bin IDs
relative to the translocation center. If "rotate" is TRUE, this combined
table is rotated relative to the left-most coordinate based on "x_bin" and
"y_bin" coordinates.

If the argument 'window_size' is provided, GRange inputs are used only to
specify the chromosome and translocation positions, then windows of the same
size are generated as appropriate. If 'window_fix' is set as "center" (the
default value), the center point of each GRange is expanded to the specified
window size (resize(fix="center")). If set to "junction", the start/end
locations of the genomic coordinates are assumed to represent the "exact"
position of the translocation, and windows are expanded in accordance with
that and with orientation. For instance, if chromosome A and B have a
translocation at 1,000,000 and 2,000,000, respectively, and both are
+ oriented, the translocation point from A and B are the End and Start
positions, respectively. If A is reversed and B is forward, these positions
are the Start and Start positions, respectively (Etc).

Retrive and combine three pixel tables describing a point
translocation: pixels from the origin sequence (side A), the destination
sequence (side B), and the pixels representing their junction (side AB).
Pixels are cached similarly to read_cooler_hdf5(), but include extra columns
"side," "x_bin,", and "y_bin." Side indicates which section of the genome
each pixel originates from, and bin columns represent arbitrary bin IDs
relative to the translocation center. If "rotate" is TRUE, this combined
table is rotated relative to the left-most coordinate based on "x_bin" and
"y_bin" coordinates.

If the argument 'window_size' is provided, GRange inputs are used only to
specify the chromosome and translocation positions, then windows of the same
size are generated as appropriate. If 'window_fix' is set as "center" (the
default value), the center point of each GRange is expanded to the specified
window size (resize(fix="center")). If set to "junction", the start/end
locations of the genomic coordinates are assumed to represent the "exact"
position of the translocation, and windows are expanded in accordance with
that and with orientation. For instance, if chromosome A and B have a
translocation at 1,000,000 and 2,000,000, respectively, and both are
+ oriented, the translocation point from A and B are the End and Start
positions, respectively. If A is reversed and B is forward, these positions
are the Start and Start positions, respectively (Etc).
}
