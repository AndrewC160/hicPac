% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_axis_title.R
\name{get_axis_title}
\alias{get_axis_title}
\title{Get axis titles}
\usage{
get_axis_title(
  pix_in,
  axis_in = "x",
  no_zero_coords = TRUE,
  append_bin_size = TRUE
)
}
\arguments{
\item{pix_in}{Matrix tibble in the format returned by read_matrix_file().}

\item{axis_in}{Which axis to return labels for. Defaults to "x", can also be "y".}

\item{no_zero_coords}{Should Zero coordinates be changed to 1? Defaults to TRUE.}
}
\description{
Given a tibble returned by read_cooler_hdf5(), axis labels.

Typical format:
"Chr1: 1,000,000-2,000,000 (1MB, 500kb bin size)" (One chromosome).
"300MB, 5MB bin size" (More than one chromosome).
}
