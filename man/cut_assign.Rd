% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_management.R
\name{cut_assign}
\alias{cut_assign}
\title{Assign cells to convert a proportion of forest.col to open invasible}
\usage{
cut_assign(pChg, ncell = NULL, assign_i = NULL, lc.df, forest.col)
}
\arguments{
\item{pChg}{Proportion of inbound cells to cut}

\item{ncell}{\code{NULL} Number of inbound grid cells. Required if cell IDs
are not supplied via \code{assign_i}}

\item{assign_i}{\code{NULL} Vector of inbound cell IDs to treat. If
\code{NULL}, then \code{ncell} IDs are sampled randomly for treatments}

\item{lc.df}{Dataframe or tibble with xy coords, land cover proportions, and
cell id info}

\item{forest.col}{Vector of forest.col column indexes within \code{lc.df}}
}
\value{
Tibble \code{id.chg} with grid and inbound indexes of cells with land
  cover change and matrix \code{mx} with the change in each forest.col
  category and the total change to be added to open habitat
}
\description{
Assign a specified number of random cells to have their forest.col habitat
cleared by a random proportion.
}
\keyword{change,}
\keyword{control,}
\keyword{cover}
\keyword{cut,}
\keyword{forest.col,}
\keyword{land}
\keyword{owners}
