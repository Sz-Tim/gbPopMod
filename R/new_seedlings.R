#' Seed germination & establishment
#'
#' This function calculates the number of new seedlings in each cell
#' @param ngrid Number of grid cells in entire map
#' @param N.seed \code{N.seed} output from \code{\link{ldd_disperse}} or
#'   \code{\link{sdd_disperse}} with grid id and number of seeds in each cell
#' @param N.sb Matrix with number of seeds in seed bank; dim=c(ngrid, tmax+1)
#' @param pr.est.ag Vector of establishment probabilities from
#'   \code{\link{cell_agg}}
#' @param pr.sb Probability of surviving a year in the seed bank
#' @param dem.st \code{Logical} denoting whether to include demographic
#'   stochasticity
#' @param bank \code{Logical} denoting whether to include a seed bank
#' @return Tibble
#' @keywords run, simulate
#' @export

new_seedlings <- function(ngrid, N.seed, N.sb, pr.est.ag, pr.sb, 
                          dem.st=F, bank=F) {
  
  require(tidyverse)
  
  N.rcrt <- rep(0, ngrid)
  if(dem.st) {
    
    N.rcrt[N.seed$id] <- rbinom(nrow(N.seed), N.seed$N, pr.est.ag[N.seed$id])
    
    if(bank) {
      N.sbEst <- rep(0, ngrid)
      
      # N_est_sb
      N.sbEst[N.seed$id] <- rbinom(nrow(N.seed), N.sb[N.seed$id], 
                                   pr.est.ag[N.seed$id])
      
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
    N.rcrt[N.seed$id] <- N.seed$N * pr.est.ag[N.seed$id,]
    
    if(bank) {
      
      # N_est_tot = N_est + N_est_sb
      N.rcrt[N.seed$id] <- (N.rcrt[N.seed$id] + 
                              N.sb[N.seed$id] * pr.est.ag[N.seed$id,]) %>% round
      
      # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
      N.sb[N.seed$id] <- ((N.sb[N.seed$id]*(1-pr.est.ag[N.seed$id,]) + 
                             N.seed$N - N.rcrt[N.seed$id]) * pr.sb) %>% round
      
    } else {
      N.sb <- rep(0, ngrid)
    }
    
  }
  return(list(N.rcrt=N.rcrt, N.sb=N.sb))
}