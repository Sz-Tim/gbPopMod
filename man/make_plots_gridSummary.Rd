% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_output.R
\name{make_plots_gridSummary}
\alias{make_plots_gridSummary}
\title{Save plots of buckthorn distribution summaries}
\usage{
make_plots_gridSummary(p, p.wd, grid.sum, cell.sum, byLC = FALSE,
  txt = NULL, w = 8, h = 6)
}
\arguments{
\item{grid.sum}{Dataframe or tibble with processed output for whole-grid
summaries}

\item{cell.sum}{Dataframe or tibble with processed output for summaries of
time-to-event by cell}

\item{byLC}{\code{FALSE} For LC-specific parameters, should separate plots be
created for each LC?}

\item{txt}{Text to append to start of plot title}

\item{w}{Figure output width (inches)}

\item{h}{Figure output height (inches)}

\item{par.wd}{Directory to store file in}
}
\value{
None
}
\description{
This function saves plots of buckthorn distribution summary statistics for
the whole grid across simulations. This includes plots for the mean and sd of
the proportion of cells occupied, of the proportion of low density cells, and
of proportion of occupied cells at carrying capacity.
}
\keyword{output}
\keyword{plots,}
\keyword{save,}
\keyword{store,}
