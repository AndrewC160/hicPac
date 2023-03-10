% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotate_points_45.R
\name{rotate_points_45}
\alias{rotate_points_45}
\title{Rotate points 45 degrees.}
\usage{
rotate_points_45(
  tb_in,
  x_grange,
  y_grange,
  x_rng = NULL,
  y_rng = NULL,
  x_col = NULL,
  y_col = NULL,
  x_col_new = NULL,
  y_col_new = NULL
)
}
\arguments{
\item{tb_in}{Tibble input.}

\item{x_grange}{GRanges object (one contiguous range) describing the limits of the x-axis.}

\item{y_grange}{GRanges object (one contiguous range) describing the limits of the y-axis.}

\item{x_rng}{Vector of length 2 with max and min values of x-range.}

\item{y_rng}{Vector of length 2 with max and min values of y-range.}

\item{x_col}{Column name of x-values.}

\item{y_col}{Column name of y-values.}

\item{x_col_new}{Name of new rotate column. Defaults to NULL, in which case x_col is replaced with the new value.}

\item{y_col_new}{Name of new rotate column. Defaults to NULL, in which case y_col is replaced with the new value.

import tibble}
}
\description{
Given a matrix with x and y columns, and a range describing an
X and y axis, rotate points in the input matrix 45 degrees and rescale them
such that the new x-axis (originally the diagonal) is the same length as the
original x-axis (originally the side).
}
