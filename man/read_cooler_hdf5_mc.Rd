% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_cooler_hdf5_mc.R
\name{read_cooler_hdf5_mc}
\alias{read_cooler_hdf5_mc}
\title{Read Cooler HDF5 (Multi-core)}
\usage{
read_cooler_hdf5_mc(
  file_cool,
  gr_range1 = NULL,
  gr_range2 = NULL,
  diag_distance = NULL,
  silent = TRUE,
  max_pixels = 6250000,
  cache_tsv = NULL,
  overwrite_cache = FALSE,
  return_table = TRUE,
  n_cores = 10,
  n_chunks = 20
)
}
\arguments{
\item{file_cool}{Name of .cool file(s) with the appropriate bin size. If multiple files are provided, functions is run recursively and tables are combined (including cache TSV).}

\item{gr_range1}{GRange of contiguous region for X dimension. Defaults to entire genome, and must be contiguous.}

\item{gr_range2}{GRange of regions for Y dimension. Defaults to values of gr_range1, and can be multiple ranges.}

\item{silent}{Should processing time messages be suppressed? Defaults to FALSE.}

\item{max_pixels}{Maximum number of pixels to retrieve without erroring out. Defaults to 6.25 million (2500x2500 grid), can also be set to Inf to disregard.}

\item{overwrite_cache}{Should cache file be overwritten? Defaults to FALSE.}

\item{return_table}{Should the table be returned? Defaults to TRUE, but if the goal is to produce the TSV cache file this option can avoid extra reading and memory requirements. Requires a cache_tsv file, if FALSE.}

\item{n_cores}{Number of cores to process data with; defaults to 10.}

\item{n_chunks}{Number of chunks to process data in; defaults to 20.}

\item{diag_distace}{Maximum distance from the diagonal (useful for filtering pixels from tables meant for rotated plots). Filter is applied AFTER fetching pixels, so does not generally help performance.}

\item{cache_table}{Filename of TSV to cache pixel results into. If provided and this file exists, this table will be read and returned. If not found, it will be created.}
}
\description{
Read .cool files as HDF5 using the R package rhdf5, and break up the job into multiple cores.If
multiple files are provided, function is run recursively on each and the results are bound into a
single table with a "cooler" column. Does not load entire cool file into memory, but instead uses
file's index to access specific rows in the pixel matrices. Includes an optional caching function
for larger data tables by providing a TSV filename for <cache_table>. Note that integer values
from Cooler may be int64, which is not natively supported by R. The rhdf5 package allows for the
use of integer64 type using the bit64 package, but as this value should only apply to genomic
coordinates counts and bin ID values are just converted to normal integers (there should not be
2 trillion counts or bins in any normal dataset). For use with data focused on the diagonal
("pyramid" plots rotated 45 degrees), the <diag_distance> argument can be supplied such that
pixels must be within this distance of the diagonal. In addition, the X axis grange (gr_range1)
is resized to be 2x the diagonal distance wider, fixed about the center. This should include all
pixels on the left and right edges which have gr_range2 values within range, but not gr_range1.
The diagonal distance filter is applied AFTER fetching pixels, so it is not an effective way to
increase retrieval speeds. It is however reflected in cached TSV files, so it will reduce file
sizes.

In some cases, it is preferable to produce the TSV cache file and NOT load its contents into
memory. In these cases, <return_table> can be set to FALSE and the pixel TSV file will be
returned; in the event this file exists it will not be read into memory.

Note that for multiple tasks to be completed at once, file locking must be disabled and thus
this function is run with the equivalent of read_cooler_hdf5(disable_file_lock=TRUE).
}
