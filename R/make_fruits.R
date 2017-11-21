#' Local fruit production
#'
#' This function calculates the number of fruits produces in ecah grid cell. It
#' assumes no fruit production before \code{age.f} and that individuals are
#' distributed among land cover types in densities relative to K.
#' @param N.t Matrix or array of abundances, with dims=c(ngrid, (lc), y.ad)
#' @param lc.mx Matrix of land cover proportions from \code{\link{cell_agg}}
#' @param fec.ag Vector of fruit per individual from \code{\link{cell_agg}}
#' @param pr.f.ag Vector of fruiting probability from \code{\link{cell_agg}}
#' @param y.ad Maximum age at maturity (i.e., \code{max(age.f)})
#' @param age.f.d \code{Logical} denoting whether \code{age.f} differs among
#'   land cover types
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @return Tibble with grid id, number of reproducing individuals, and total
#'   fruit produces within the cell
#' @keywords fruit, reproduction, fecundity
#' @export

make_fruits <- function(N.t, lc.mx, fec.ag, pr.f.ag, y.ad, age.f.d, dem.st=F) {
  
  require(tidyverse)

  # calculate N.mature in each LC in each cell
  if(age.f.d) {
    N.mature <- rowSums(N.t[,,y.ad])
  } else {
    N.mature <- N.t[,y.ad]
  }
  if(dem.st) {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.rpr = rbinom(n(), N.mature[id],
                            prob=pr.f.ag[id]),
             N.fruit = rpois(n(), 
                             lambda=N.rpr*fec.ag[id])) %>% 
      filter(N.fruit > 0)
  } else {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.rpr=(N.mature[id]) * pr.f.ag[id,],
             N.fruit=(N.rpr * fec.ag[id,]) %>% round) %>% 
      filter(N.fruit > 0)
  }
  return(N.f)
}