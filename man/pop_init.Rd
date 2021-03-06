% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_aux.R
\name{pop_init}
\alias{pop_init}
\title{Initialize populations randomly}
\usage{
pop_init(ngrid, g.p, lc.df, p.0 = NULL, N.0 = NULL)
}
\arguments{
\item{ngrid}{Number of grid cells in entire map}

\item{g.p}{Named list of global parameters}

\item{lc.df}{Dataframe or tibble with xy coords, land cover proportions, and
cell id info}

\item{p.0}{\code{NULL} Cell IDs with populations at t=0}

\item{N.0}{\code{NULL} Initial adult abundance in cell}
}
\value{
Matrix or array of initial abundances with \code{dim=c(ngrid, (n.lc),
  m.max)}
}
\description{
Initialize populations randomly with populated cells containing adults at 50%
K and juveniles at 10% K
}
\keyword{initialize,}
\keyword{set}
\keyword{up}
