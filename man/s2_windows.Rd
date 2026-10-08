% Generated documentation for R/s2_windows.R
\name{s2_windows}
\alias{s2_windows}
\docType{data}
\title{Sentinel-2 spectral-index windows}
\format{
A named list of three numeric matrices, each with 45 rows and 40 columns.
Rows are grouped by index in the order NBR, MSI, CRSWIR, CRRE, NMDI,
with nine cells per index. The fifth cell in each group is the focal cell
(rows 5, 14, 23, 32 and 41). Columns are ordered by year and within-year
position (1 through 4), with time identifiers 201601, 201602, ..., 202504.
The list names retain the source-window identifiers: window_01, window_02
and window_04.
}
\source{
Input matrices supplied by Donato Morresi, extracted from a FORCE data cube.
}
\usage{
s2_windows
}
\description{
Three input matrices extracted from FORCE data using a 3 by 3 spatial kernel.
Each matrix contains five spectral indices observed four times per year
from 2016 through 2025, corresponding to June through September.
Missing observations are retained for gap filling by \code{\link[=hilandyn2_win]{hilandyn2_win()}}.
}
\seealso{
\code{\link[=hilandyn2_win]{hilandyn2_win()}} for an executable comparison
with and without seasonal adjustment.
}
\keyword{datasets}
