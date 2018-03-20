#' Set short distance dispersal probabilities
#'
#' Assign base dispersal probabilities from each cell. Each layer \code{k} in
#' \code{[1:i, 1:j, k, source.id]} is the SDD neighborhood for cell n. k=1
#' contains pr(SDD | source.id,i,j); k=2 contains the grid id for each cell in
#' the neighborhood.
#' @param ncell Number of inbound grid cells
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param g.p Named list of global parameters
#' @param lc.new Vector of indexes of cells that need to have their SDD
#'   neighborhood recalculated; defaults to \code{NULL} and calculates SDD
#'   neighborhoods for all cells
#' @param edges \code{"wall"} Boundary behavior
#' @param lc.col \code{4:9} Column indexes for land cover proportions
#' @return Array with dim(disp.rows, disp.cols, 2, ncell) where the third
#'   dimension contains grid id's for the neighborhood or probabilities to each
#'   target cell
#' @keywords sdd, dispersal, probability, probabilities
#' @export

sdd_set_probs <- function(ncell, lc.df, g.p, lc.new=NULL, 
                          edges="wall", lc.col=4:9) {
  
  library(purrr); library(tidyverse); library(fastmatch)
  
  # unpack parameters
  sdd.max <- g.p$sdd.max
  sdd.rate <- g.p$sdd.rate
  bird.hab <- g.p$bird.hab
  
  # initialize landscape & storage objects
  n.x <- range(lc.df[,1])
  n.y <- range(lc.df[,2])
  nbr <- 2 * sdd.max + 1
  sdd.i <- array(0, dim=c(nbr, nbr, 2, ncell))
  bird.hab.E <- as.matrix(lc.df[,lc.col]) %*% (bird.hab %>% divide_by(sum(.)))
  if(edges=="wall") bird.hab.E[!lc.df$inbd] <- 0
  
  # generate default dispersal probability matrix
  d.pr <- matrix(0, nbr, nbr)
  ctr <- sdd.max + 1  # center index (i=j) for square mx
  for(i in 1:nbr) {
    for(j in i:nbr) {
      d.pr[i,j] <- dexp((i-ctr)^2 + (j-ctr)^2 - 0.5, sdd.rate)
      if( sqrt((i-ctr)^2 + (j-ctr)^2) > sdd.max ) {
        d.pr[i,j] <- 0
      }
      d.pr[j,i] <- d.pr[i,j]
    }
  }
  d.pr <- d.pr/sum(d.pr)
  
  # pair cell IDs for each neighborhood; indexes match neighborhood matrix
  if(is.null(lc.new)) {
    if(edges=="none") {
      xx <- map(lc.df$x, ~seq(.-sdd.max, .+sdd.max))
      yy <- map(lc.df$y, ~seq(.-sdd.max, .+sdd.max))
    } else {
      xx <- map(lc.df$x[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
      yy <- map(lc.df$y[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
    }
  } else {
    xx <- map(lc.df$x[lc.df$id %in% lc.new$id], ~seq(.-sdd.max, .+sdd.max))
    yy <- map(lc.df$y[lc.df$id %in% lc.new$id], ~seq(.-sdd.max, .+sdd.max))
  }
  
  # create lists of on-map xy neighborhood ranges
  n.ix <- map(xx, ~.[.>=n.x[1] & .<=n.x[2]])
  n.iy <- map(yy, ~.[.>=n.y[1] & .<=n.y[2]])
  
  # generate all xy combinations & neighborhood matrix indices
  if(is.null(lc.new)) cat("  generating neighborhoods...\n")
  n.i <- map2(n.ix, n.iy, expand_v) 
  n.x <- map2(xx, n.ix, `%fin%`) %>% map(which) %>% map(range)
  n.y <- map2(yy, n.iy, `%fin%`) %>% map(which) %>% map(range)
  
  # match xy combinations with cell IDs
  if(is.null(lc.new)) cat("  determining neighborhood cell IDs...\n")
  c.i <- map(n.i, ~fmatch(., lc.df$x_y))

  if(is.null(lc.new)) cat("  calculating probabilities...\n")
  for(n in 1:ncell) {
    # find cell ID for each cell in neighborhood
    sdd.i[n.y[[n]][1]:n.y[[n]][2],
          n.x[[n]][1]:n.x[[n]][2],2,n] <- matrix(c.i[[n]], 
                                                 ncol=diff(n.x[[n]])+1,
                                                 byrow=TRUE)
    # weight by bird habitat preference & set cell ID to 0 if pr(target) == 0 
    ib <- sdd.i[,,2,n] != 0
    sdd.i[,,1,n][ib] <- d.pr[ib] * bird.hab.E[sdd.i[,,2,n][ib]]
    sdd.i[,,2,n][sdd.i[,,1,n]==0] <- 0
    
    # progress update
    if(n %% 5000 == 0) {
      if(is.null(lc.new)) cat("  finished cell", n, "\n")
    }
  }
  if(is.null(lc.new)) cat("  finished:", n, "cells\n")
  sdd.i[,,1,] <- apply(sdd.i[,,1,], 3, function(x) x/sum(x))
  return(sdd.i)
}




#' Short distance dispersal
#'
#' This function calculates the number of viable seeds deposited within the
#' short distance dispersal neighborhood of each cell, accounting for the
#' distance from the source cell, bird habitat preference, the proportion of
#' fruits eaten by birds, and seed viability post-digestion.
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.in} indexes only inbound cells
#' @param N.f Tibble of fruits produced in each cell output from
#'   \code{\link{make_fruits}}
#' @param pr.eat.E Vector of proportion of fruits eaten by birds from
#'   \code{\link{cell_agg}}
#' @param pr.s.bird Proportion of viable seeds post-digestion
#' @param sdd.pr Array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell) output
#'   from \code{\link{sdd_set_probs}}
#' @param sdd.rate Rate parameter for SDD exponential kernel
#' @param sdd.st \code{Logical} denoting whether to implement short distance
#'   dispersal stochastically
#' @param edges Character taking the value of one of: \code{"wall", "sink",
#'   "none"} where \code{"wall"} results in a dispersal probability of 0 for all
#'   out-of-bound cells with no populations modeled, \code{"sink"} results in
#'   dispersal of seeds to out-of-bound cells but no populations modeled, and
#'   \code{"none"} results in dispersal of seeds and populations modeled
#' @return Tibble with grid id and number of seeds in each cell
#' @keywords dispersal, SDD
#' @export

sdd_disperse <- function(id.i, N.f, pr.eat.E, pr.s.bird, 
                         sdd.pr, sdd.rate, sdd.st=F, edges="wall") {
  
  library(tidyverse); library(magrittr)
  
  # calculate seeds deposited within source cell vs emigrants
  N.source <- N.f %>%
    mutate(N.produced=(2.3*N.fruit),
           N.emig=N.produced*(1-pexp(.5,sdd.rate))*pr.eat.E[id,],
           N.drop=N.produced-N.emig) %>%
    mutate(N.emig=N.emig*pr.s.bird,
           id.in=id.i$id.in[id])
  N.seed <- N.source %>% select(id, N.drop) %>% rename(N.dep=N.drop)
  
  if(sdd.st) {
    N.seed$N.dep <- round(N.seed$N.dep)
    if(edges=="none") {
      SDD.sd <- unlist(apply(N.source, 1,
                             function(x) sample(sdd.pr[,,2,x[1]], x[5], 
                                                replace=TRUE,
                                                prob=sdd.pr[,,1,x[1]])))
    } else {
      SDD.sd <- unlist(apply(N.source, 1,
                             function(x) sample(sdd.pr[,,2,x[7]], x[5], 
                                                replace=TRUE,
                                                prob=sdd.pr[,,1,x[7]])))
    }
    SDD.dep <- tabulate(SDD.sd)  # vector of counts for 1:max(SDD.sd)
    SDD.nonzero <- SDD.dep > 0  # cell id's with N.dep > 0
    N.seed <- add_row(N.seed, 
                      id=which(SDD.nonzero), 
                      N.dep=SDD.dep[SDD.nonzero])
  } else {
    # assign emigrants to target cells & sum within each cell
    if(edges=="none") {
      N.seed %<>% 
        add_row(id=apply(N.source, 1, 
                         function(x) c(sdd.pr[,,2,x[1]])) %>% c, 
                N.dep=apply(N.source, 1, 
                            function(x) c(x[5] * sdd.pr[,,1,x[1]])) %>% c) %>%
        filter(N.dep > 0)
    } else {
      N.seed %<>% 
        add_row(id=apply(N.source, 1, 
                         function(x) c(sdd.pr[,,2,x[7]])) %>% c, 
                N.dep=apply(N.source, 1, 
                            function(x) c(x[5] * sdd.pr[,,1,x[7]])) %>% c) %>%
        filter(N.dep > 0)
    }
  }
  
  N.seed %<>%
    group_by(id) %>% 
    summarise(N=sum(N.dep) %>% round) %>%
    filter(N > 0)
  if(edges=="wall") N.seed %<>% filter(!is.na(id.i$id.in[id]))
  
  return(N.seed)
}




#' Long distance dispersal
#'
#' This function assigns n.ldd random long distance dispersal events across the
#' landcape. A single established seed is added to each target cell
#' @param ncell Number of inbound grid cells
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.in} indexes only inbound cells
#' @param N.rcrt \code{N.rcrt} output from \code{\link{new_seedlings}} with grid
#'   id and number of new recruits in each cell
#' @param n.ldd Number of long distance dispersal events per time step
#' @return Tibble with grid id and number of seeds
#' @keywords LDD, dispersal
#' @export

ldd_disperse <- function(ncell, id.i, N.rcrt, n.ldd) {
  
  ldd.id <- id.i$id[which(id.i$id.in %in% sample(1:ncell, n.ldd, replace=T))]
  N.rcrt[ldd.id] <- N.rcrt[ldd.id] + 1

  return(N.rcrt)
}







#' Short distance dispersal: lambda-based
#'
#' This function calculates the number of immigrants to each cell, accounting
#' for distance from source cell and bird habitat preferences. It uses the
#' output from simple lambda-based population growth rather than the demographic
#' version. The output values account for carrying capacity.
#' @param ncell Number of inbound grid cells
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.in} indexes only inbound cells
#' @param sdd.pr Array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell) output
#'   from \code{\link{sdd_set_probs}}
#' @param sdd.rate Rate parameter for SDD exponential kernel
#' @param K.E Carrying capacity of each cell output from \code{\link{cell_E}}
#' @param sdd.st \code{Logical} denoting whether to implement short distance
#'   dispersal stochastically
#' @return Tibble with grid id and number of individuals limited by K
#' @keywords SDD, dispersal, lambda
#' @export

sdd_lambda <- function(N.new, id.i, sdd.pr, sdd.rate, K.E, sdd.st=F) {
  # Calculate (N.arrivals | N.new, sdd.probs)
  # Accounts for distance from source cell & bird habitat preference
  # Returns dataframe with total population sizes.
  
  N.source <- N.new %>% filter(N.new > 0) %>% mutate(id.in=id.i$id.in[id])
  N.emig <- tibble(id=id.i$id, 
                   N=N.new$N.pop.upd[match(id.i$id, N.new$id)]) %>% 
    mutate(id.in=id.i$id.in[id])
  if(sdd.st) {
    SDD.sd <- unlist(apply(N.source, 1,
                           function(x) sample(sdd.pr[,,2,x[5]], 
                                              x[3] * (1-pexp(.5,sdd.rate)), 
                                              replace=TRUE,
                                              prob=sdd.pr[,,1,x[5]])))
    SDD.dep <- tabulate(SDD.sd)  # vector of counts for 1:max(SDD.sd)
    SDD.nonzero <- SDD.dep > 0  # cell id's with N.dep > 0
    N.emig %<>% add_row(id=which(SDD.nonzero), 
                        N=SDD.dep[SDD.nonzero],
                        id.in=id.i$id.in[match(id, id.i$id)])
  } else {
    N.emig %<>%
      add_row(id=apply(N.source, 1, function(x) c(sdd.pr[,,2,x[5]])) %>% c, 
              N=apply(N.source, 1, function(x) c(x[3] * (1-pexp(.5,sdd.rate)) * 
                                                   sdd.pr[,,1,x[5]])) %>% c,
              id.in=id.i$id.in[match(id, id.i$id)]) 
  }
  N.emig %<>% filter(id != 0) %>% group_by(id) %>%
    summarise(N=sum(N, na.rm=T))
  N.emig$N <- round(pmin(K.E[N.emig$id,], N.emig$N))
  
  return(N.emig)
}

