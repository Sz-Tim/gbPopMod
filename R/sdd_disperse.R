#' Short distance dispersal
#'
#' This function calculates the number of viable seeds deposited within the
#' short distance dispersal neighborhood of each cell, accounting for the
#' distance from the source cell, bird habitat preference, the proportion of
#' fruits eaten by birds, and seed viability post-digestion.
#' @param id_i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id_inbd} indexes only inbound cells
#' @param N.f Tibble of fruits produced in each cell output from
#'   \code{\link{make_fruits}}
#' @param pr.eat.ag Vector of proportion of fruits eaten by birds from
#'   \code{\link{cell_agg}}
#' @param pr.s.bird Proportion of viable seeds post-digestion
#' @param sdd.pr Array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell) output
#'   from \code{\link{sdd_set_probs}}
#' @param sdd.rate Rate parameter for SDD exponential kernel
#' @param sdd.st \code{Logical} denoting whether to implement short distance
#'   dispersal stochastically
#' @return Tibble with grid id and number of seeds in each cell
#' @keywords dispersal, SDD
#' @export

sdd_disperse <- function(id_i, N.f, pr.eat.ag, pr.s.bird, 
                   sdd.pr, sdd.rate, sdd.st=F) {
  
  require(tidyverse); require(magrittr)
  
  # calculate seeds deposited within source cell vs emigrants
  N.source <- N.f %>%
    mutate(N.produced=(2.3*N.fruit),
           N.emig=N.produced*(1-pexp(.5,sdd.rate))*pr.eat.ag[id,],
           N.drop=N.produced-N.emig) %>%
    mutate(N.emig=N.emig*pr.s.bird,
           id_inbd=id_i$id_inbd[id])
  N.seed <- N.source %>% select(id, N.drop) %>% rename(N.dep=N.drop)
  
  if(sdd.st) {
    N.seed$N.dep <- round(N.seed$N.dep)
    SDD_sd <- unlist(apply(N.source, 1,
                           function(x) sample(sdd.pr[,,2,x[7]], x[5], 
                                              replace=TRUE,
                                              prob=sdd.pr[,,1,x[7]])))
    SDD_dep <- tabulate(SDD_sd)  # vector of counts for 1:max(SDD_sd)
    SDD_nonzero <- SDD_dep > 0  # cell id's with N_dep > 0
    N.seed <- add_row(N.seed, 
                      id=which(SDD_nonzero), 
                      N.dep=SDD_dep[SDD_nonzero])
  } else {
    # assign emigrants to target cells & sum within each cell
    N.seed %<>% 
      add_row(id=apply(N.source, 1, 
                       function(x) c(sdd.pr[,,2,x[7]])) %>% c, 
              N.dep=apply(N.source, 1, 
                          function(x) c(x[5] * sdd.pr[,,1,x[7]])) %>% c) %>%
      filter(N.dep > 0)
  }
  
  N.seed %<>%
    group_by(id) %>% 
    summarise(N=sum(N.dep) %>% round) %>%
    filter(!is.na(id_i$id_inbd[id]) & N > 0)
  
  return(N.seed)
}