#' Local fruit production
#'
#' This function calculates the number of fruits produces in ecah grid cell. It
#' assumes no fruit production before \code{age.f} and that individuals are
#' distributed among land cover types in densities relative to K.
#' @param N.t Matrix or array of abundances, with dims=c(ngrid, (lc), y.ad)
#' @param lc.mx Matrix of land cover proportions from \code{\link{cell_E}}
#' @param fec.E Vector of fruit per individual from \code{\link{cell_E}}
#' @param pr.f.E Vector of fruiting probability from \code{\link{cell_E}}
#' @param y.ad Maximum age at maturity (i.e., \code{max(age.f)})
#' @param age.f.d \code{Logical} denoting whether \code{age.f} differs among
#'   land cover types
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @return Tibble with grid id, number of reproducing individuals, and total
#'   fruit produces within the cell
#' @keywords fruit, reproduction, fecundity
#' @export

make_fruits <- function(N.t, lc.mx, fec.E, pr.f.E, y.ad, age.f.d, dem.st=F) {
  
  library(tidyverse)
  
  # calculate N.mature in each LC in each cell
  if(age.f.d) {
    N.mature <- rowSums(N.t[,,y.ad])
  } else {
    N.mature <- N.t[,y.ad]
  }
  if(dem.st) {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.rpr = rbinom(n(), N.mature[id], prob=pr.f.E[id]),
             N.fruit = rpois(n(), lambda=N.rpr*fec.E[id])) %>% 
      filter(N.fruit > 0)
  } else {
    N.f <- tibble(id = which(N.mature>0)) %>%
      mutate(N.rpr=(N.mature[id]) * pr.f.E[id,],
             N.fruit=(N.rpr * fec.E[id,]) %>% round) %>% 
      filter(N.fruit > 0)
  }
  return(N.f)
}




#' Seed germination & establishment
#'
#' This function calculates the number of new seedlings in each cell
#' @param ngrid Number of grid cells in entire map
#' @param N.seed \code{N.seed} output from \code{\link{ldd_disperse}} or
#'   \code{\link{sdd_disperse}} with grid id and number of seeds in each cell
#' @param N.sb Matrix with number of seeds in seed bank; dim=c(ngrid, tmax+1)
#' @param pr.est.E Vector of establishment probabilities from
#'   \code{\link{cell_E}}
#' @param pr.sb Probability of surviving a year in the seed bank
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @param bank \code{Logical} denoting whether to include a seed bank
#' @return Tibble
#' @keywords run, simulate
#' @export

new_seedlings <- function(ngrid, N.seed, N.sb, pr.est.E, pr.sb, 
                          dem.st=F, bank=F) {
  
  library(tidyverse)
  
  N.rcrt <- rep(0, ngrid)
  if(dem.st) {
    
    N.rcrt[N.seed$id] <- rbinom(nrow(N.seed), N.seed$N, pr.est.E[N.seed$id])
    
    if(bank) {
      N.sbEst <- rep(0, ngrid)
      
      # N_est_sb
      N.sbEst[N.seed$id] <- rbinom(nrow(N.seed), N.sb[N.seed$id], 
                                   pr.est.E[N.seed$id])
      # N_est_tot = N_est + N_est_sb
      N.rcrt[N.seed$id] <- N.rcrt[N.seed$id] + N.sbEst[N.seed$id]
      # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
      N.sb[N.seed$id] <- rbinom(nrow(N.seed),
                                N.sb[N.seed$id] + N.seed$N - N.rcrt[N.seed$id],
                                pr.sb)
    } else {
      N.sb <- rep(0, ngrid)
    }
    
  } else {
    
    # N_est = N_seed * p(est)
    N.rcrt[N.seed$id] <- N.seed$N * pr.est.E[N.seed$id,]
    
    if(bank) {
      
      # N_est_tot = N_est + N_est_sb
      N.rcrt[N.seed$id] <- (N.rcrt[N.seed$id] + 
                              N.sb[N.seed$id] * pr.est.E[N.seed$id,]) %>% round
      # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
      N.sb[N.seed$id] <- ((N.sb[N.seed$id]*(1-pr.est.E[N.seed$id,]) + 
                             N.seed$N - N.rcrt[N.seed$id]) * pr.sb) %>% round
    } else {
      N.sb <- rep(0, ngrid)
    }
    
  }
  return(list(N.rcrt=N.rcrt, N.sb=N.sb))
}





#' Local population growth: lambda-based (a la Merow 2011)
#'
#' This function calculates change in population size in each cell with
#' lambda-based proportional growth.
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
