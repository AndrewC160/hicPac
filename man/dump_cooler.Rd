% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dump_cooler.R
\name{dump_cooler}
\alias{dump_cooler}
\title{Dump cooler}
\usage{
dump_cooler(
  gr_range1 = NULL,
  gr_range2 = NULL,
  file_cool,
  silent = FALSE,
  exe_cooler = "/N/u/aclugston/anaconda3/envs/EagleC/bin/cooler"
)
}
\arguments{
\item{gr_range1}{GRange of contiguous region for X dimension. Defaults to entire genome.}

\item{gr_range2}{GRange of contiguous region for Y dimension. Defaults to values of gr_range1.}

\item{file_cool}{Filename of .cool file with the appropriate bin size.}

\item{exe_cooler}{Filename of cooler executable. Defaults to '~/anaconda3/envs/EagleC/bin/cooler'.}
}
\description{
Dump cooler contents from a given region, either in one square (same range
in two dimensions) or in an asymmetric rectangle (one range for X dimension,
another for Y).
}
