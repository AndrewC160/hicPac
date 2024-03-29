% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlf_parse_assembly.R
\name{nlf_parse_assembly}
\alias{nlf_parse_assembly}
\title{NeoLoopFinder parse assembly}
\usage{
nlf_parse_assembly(txt_in)
}
\arguments{
\item{txt_in}{Text line from assembly text file.}
}
\description{
Given a line of text from a NeoLoopFinder assembly file, parse it into a table
which is organized as one line per segment with the assembly ID listed in an 
"assembly" column. Strand orientation is represented in the "strand" column, 
and each segment is labeled as "seg_1" etc. Columns "type_start" and "type_end"
indicate which type of SV produced each junction in the assembly (start on the
left and end on the right), and in the case of the left- or right-most ends of
the assembly these will be labeled as "left_bound" or "right_bound",
respectively.
}
