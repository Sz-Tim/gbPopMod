% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_dispersal.R
\name{sdd_disperse}
\alias{sdd_disperse}
\title{Short distance dispersal}
\usage{
sdd_disperse(id.i, Fr, gamma, p.c.E, s.c, sdd.sp, sdd.rate, sdd.st = T,
  edges = "wall")
}
\arguments{
\item{id.i}{Tibble matching cell IDs. \code{id} indexes on the entire grid
while \code{id.in} indexes only inbound cells}

\item{Fr}{Tibble of fruits produced in each cell output from
\code{\link{make_fruits}}}

\item{gamma}{Scalar: mean number of seeds per fruit}

\item{p.c.E}{Vector of proportion of fruits eaten by birds from
\code{\link{cell_E}}}

\item{s.c}{Proportion of viable seeds post-digestion}

\item{sdd.sp}{List of named vectors with SDD probabilities for each cell
output from \code{\link{sdd_set_probs}}}

\item{sdd.rate}{Rate parameter for SDD exponential kernel}

\item{sdd.st}{\code{Logical} denoting whether to implement short distance
dispersal stochastically}

\item{edges}{Character taking the value of one of: \code{"wall", "sink",
"none"} where \code{"wall"} results in a dispersal probability of 0 for all
out-of-bound cells with no populations modeled, \code{"sink"} results in
dispersal of seeds to out-of-bound cells but no populations modeled, and
\code{"none"} results in dispersal of seeds and populations modeled}
}
\value{
Tibble with grid id and number of seeds in each cell
}
\description{
Calculate the number of viable seeds deposited within the short distance
dispersal neighborhood of each cell, accounting for the distance from the
source cell, bird habitat preference, the proportion of fruits eaten by birds,
and seed viability post-digestion.
}
\keyword{SDD}
\keyword{dispersal,}
