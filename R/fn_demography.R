#' Local fruit production
#'
#' Calculate the number of fruits produces in ecah grid cell. It assumes no
#' fruit production before \code{m} and that individuals are distributed among
#' land cover types in densities relative to K.
#' @param N.t Matrix or array of abundances, with dims=c(ngrid, (lc), m.max)
#' @param lc.mx Matrix of land cover proportions from \code{\link{cell_E}}
#' @param mu.E Vector of fruit per individual from \code{\link{cell_E}}
#' @param p.f.E Vector of fruiting probability from \code{\link{cell_E}}
#' @param m.max Maximum age at maturity (i.e., \code{max(m)})
#' @param m.d \code{Logical} denoting whether \code{m} differs among land cover
#'   types
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @return Tibble with grid id, number of reproducing individuals, and total
#'   fruit produces within the cell
#' @keywords fruit, reproduction, fecundity
#' @export

make_fruits <- function(N.t, lc.mx, mu.E, p.f.E, m.max, m.d, dem.st=F) {
  
  library(tidyverse)
  
  # calculate N.mature in each LC in each cell
  if(m.d) {
    N.mature <- rowSums(N.t[,,m.max])
  } else {
    N.mature <- N.t[,m.max]
  }
  if(dem.st) {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.ad=N.mature[id],
             N.rpr = rbinom(n(), N.mature[id], prob=p.f.E[id]),
             N.fruit = rpois(n(), lambda=N.rpr*mu.E[id])) %>% 
      filter(N.fruit > 0)
  } else {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.ad=N.mature[id],
             N.rpr=(N.mature[id] * p.f.E[id,]) %>% round,
             N.fruit=(N.rpr * mu.E[id,]) %>% round) %>% 
      filter(N.fruit > 0)
  }
  return(N.f)
}




#' Seed germination & establishment
#'
#' Calculate the number of new seedlings in each cell
#' @param ngrid Number of grid cells in entire map
#' @param N.seed \code{N.seed} output from \code{\link{ldd_disperse}} or
#'   \code{\link{sdd_disperse}} with grid id and number of seeds in each cell
#' @param B Matrix with number of seeds in seed bank; dim=c(ngrid, tmax+1)
#' @param p.E Vector of establishment probabilities from
#'   \code{\link{cell_E}}
#' @param g.D Probability of germinating in same year as produced
#' @param g.B Probability of germinating from seed bank
#' @param s.B Probability of surviving a year in the seed bank
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @param bank \code{Logical} denoting whether to include a seed bank
#' @return Tibble
#' @keywords run, simulate
#' @export

new_seedlings <- function(ngrid, N.seed, B, p.E, g.D, g.B, s.B, 
                          dem.st=F, bank=F) {
  
  library(tidyverse)
  
  M.D <- germ.D <- rep(0, ngrid)
  if(dem.st) {
    M.0[N.seed$id] <- rbinom(nrow(N.seed), N.seed$N, p.E[N.seed$id])
    if(bank) {
      N.sbEst <- rep(0, ngrid)
      # N_est_sb
      N.sbEst[N.seed$id] <- rbinom(nrow(N.seed), B[N.seed$id], 
                                   p.E[N.seed$id])
      # N_est_tot = N_est + N_est_sb
      M.0[N.seed$id] <- M.0[N.seed$id] + N.sbEst[N.seed$id]
      # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
      B[N.seed$id] <- rbinom(nrow(N.seed),
                                B[N.seed$id] + N.seed$N - M.0[N.seed$id],
                                s.B)
    } else {
      B <- rep(0, ngrid)
    }
    
  } else {
    # N_est = N_seed * g.D * p(est)
    germ.D[N.seed$id] <- round(N.seed$N * g.D)
    M.D <- round(germ.D * p.E)
    if(bank) {
      # N_est_sb = B * g.B * p(est)
      germ.B <- round(B * g.B)
      M.B <- round(germ.B * p.E)
      # N_est_tot = N_est + N_est_sb
      M.0 <- M.D + M.B
      # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
      B <- B - germ.B
      B[N.seed$id] <- B[N.seed$id] + N.seed$N - germ.D[N.seed$id]
      B <- round(B * s.B)
    } else {
      B <- rep(0, ngrid)
      M.0 <- M.D
    }
  }
  return(list(M.0=M.0, B=B))
}





#' Local population growth: lambda-based (a la Merow 2011)
#'
#' Calculate change in population size in each cell with lambda-based
#' proportional growth.
#' @param N.t Vector of abundances with \code{length=ngrid}
#' @param lambda.E Vector with lambda for each cell with \code{length=ngrid}
#' @param sdd.rate Rate parameter for SDD exponential kernel
#' @return Sparse tibble with grid id, starting population size, change in
#'   population size, and updated population size accounting for emigrants
#' @keywords lambda, growth
#' @export

grow_lambda <- function(N.t, lambda.E, sdd.rate) {
  
  N.id <- which(N.t>0)
  lam.id <- lambda.E[N.id]
  N.new <- tibble(id = N.id) %>%
    mutate(N.pop=N.t[N.id],
           N.new=N.pop * (lam.id-1),
           N.pop.upd=N.pop + (lam.id>=1) * N.new * pexp(0.5, sdd.rate) +
             (lam.id<1) * N.new)
  
  return(N.new)
}
