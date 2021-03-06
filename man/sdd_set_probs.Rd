% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_dispersal.R
\name{sdd_set_probs}
\alias{sdd_set_probs}
\title{Set short distance dispersal probabilities}
\usage{
sdd_set_probs(ncell, lc.df, g.p, edges = "wall", lc.col = 4:9,
  verbose = F)
}
\arguments{
\item{ncell}{Number of inbound grid cells}

\item{lc.df}{Dataframe or tibble with xy coords, land cover proportions, and
cell id info}

\item{g.p}{Named list of global parameters}

\item{edges}{\code{"wall"} Boundary behavior}

\item{lc.col}{\code{4:9} Column indexes for land cover proportions}

\item{verbose}{\code{FALSE} Give updates for number of cells completed?}
}
\value{
List with full neighborhoods,i, sparse representation, sp, and sparse
  dataframe sp.df. The full array has dim(disp.rows, disp.cols, 2, ncell)
  where the third dimension contains grid id's for the neighborhood or
  probabilities to each target cell. The sparse representation contains a
  list with containing the cells dispersing into each cell and a list with
  the associated probabilities. The sparse dataframe is a dataframe with a
  row for each non-zero i-j pair with columns for i, j, and dispersal
  probability
}
\description{
Assign base dispersal probabilities from each cell. Each layer \code{k} in
\code{[1:i, 1:j, k, source.id]} is the SDD neighborhood for cell n. k=1
contains pr(SDD | source.id,i,j); k=2 contains the grid id for each cell in
the neighborhood.
}
\keyword{dispersal,}
\keyword{probabilities}
\keyword{probability,}
\keyword{sdd,}
