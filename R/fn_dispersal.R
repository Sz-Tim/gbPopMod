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
#' @param edges \code{"wall"} Boundary behavior
#' @param lc.col \code{4:9} Column indexes for land cover proportions
#' @return List with full neighborhoods,i, sparse representation, sp, and sparse
#'   dataframe sp.df. The full array has dim(disp.rows, disp.cols, 2, ncell)
#'   where the third dimension contains grid id's for the neighborhood or
#'   probabilities to each target cell. The sparse representation contains a
#'   list with containing the cells dispersing into each cell and a list with
#'   the associated probabilities. The sparse dataframe is a dataframe with a
#'   row for each non-zero i-j pair with columns for i, j, and dispersal
#'   probability
#' @param verbose \code{FALSE} Give updates for number of cells completed?
#' @keywords sdd, dispersal, probability, probabilities
#' @export

sdd_set_probs <- function(ncell, lc.df, g.p, 
                          edges="wall", lc.col=4:9, verbose=F) {
  
  library(tidyverse); library(magrittr); library(fastmatch)
  
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
  if(edges=="none") {
    xx <- map(lc.df$x, ~seq(.-sdd.max, .+sdd.max))
    yy <- map(lc.df$y, ~seq(.-sdd.max, .+sdd.max))
  } else {
    xx <- map(lc.df$x[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
    yy <- map(lc.df$y[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
  }
  
  # create lists of on-map xy neighborhood ranges
  n.ix <- map(xx, ~.[.>=n.x[1] & .<=n.x[2]])
  n.iy <- map(yy, ~.[.>=n.y[1] & .<=n.y[2]])
  
  # generate all xy combinations & neighborhood matrix indices
  if(verbose) cat("  generating neighborhoods...\n")
  n.i <- map2(n.ix, n.iy, expand_v) 
  n.x <- map2(xx, n.ix, `%fin%`) %>% map(which) %>% map(range)
  n.y <- map2(yy, n.iy, `%fin%`) %>% map(which) %>% map(range)
  
  # match xy combinations with cell IDs
  if(verbose) cat("  determining neighborhood cell IDs...\n")
  c.i <- map(n.i, ~fmatch(., lc.df$x_y))


  if(verbose) cat("  calculating probabilities...\n")
  if(verbose) pb <- txtProgressBar(min=1, max=ncell, style=3)

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
    

    if(verbose) setTxtProgressBar(pb, n)
  }
  if(verbose) {cat("  finished:", n, "cells\n"); close(pb)}

  sdd.i[,,1,] <- apply(sdd.i[,,1,], 3, function(x) x/sum(x))
  sdd.sparse.ls <- vector("list", ncell)
  for(n in 1:ncell) {
    sdd.sparse.ls[[n]] <- c(sdd.i[,,1,n])
    names(sdd.sparse.ls[[n]]) <- c(sdd.i[,,2,n])
    sdd.sparse.ls[[n]] <- sdd.sparse.ls[[n]][sdd.sparse.ls[[n]]>0]
  }
  sdd.sparse <- data.frame(i.idin=rep(1:ncell, 
                                 times=map_int(sdd.sparse.ls, length)),
                           j.id=unlist(map(sdd.sparse.ls, ~as.numeric(names(.)))),
                           pr=unlist(sdd.sparse.ls))
  sdd.sparse$j.idin <- lc.df$id.in[match(sdd.sparse$j.id, lc.df$id)]
  return(list(i=sdd.i, sp=sdd.sparse.ls, sp.df=sdd.sparse))
}




#' Update short distance dispersal probabilities
#'
#' Update dispersal probabilities for each affected cell following land cover
#' changes. Because SDD probabilities are weighted by land cover based on bird
#' habitat preferences, the probabilities must be recalculated for each cell
#' that contains an affected cell as its target. Each layer \code{k} in
#' \code{[1:i, 1:j, k, source.id]} is the SDD neighborhood for cell n. k=1
#' contains pr(SDD | source.id,i,j); k=2 contains the grid id for each cell in
#' the neighborhood.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param g.p Named list of global parameters
#' @param sdd.alter Vector of indexes of cells that need to have their SDD
#'   neighborhood recalculated
#' @param sdd.0 Original (non-sparse) SDD neighborhoods generated by
#'   \link{sdd_set_probs}; i.e., sdd_set_probs()$i
#' @param edges \code{"wall"} Boundary behavior
#' @param lc.col \code{4:9} Column indexes for land cover proportions
#' @return List with full neighborhoods,i, sparse representation, sp, and sparse
#'   dataframe sp.df. The full array has dim(disp.rows, disp.cols, 2, ncell)
#'   where the third dimension contains grid id's for the neighborhood or
#'   probabilities to each target cell. The sparse representation contains a
#'   list with containing the cells dispersing into each cell and a list with
#'   the associated probabilities. The sparse dataframe is a dataframe with a
#'   row for each non-zero i-j pair with columns for i, j, and dispersal
#'   probability
#' @export

sdd_update_probs <- function(lc.df, g.p, sdd.alter, sdd.0, lc.col=4:9) {
  
  library(tidyverse); library(magrittr)
  
  # unpack parameters
  sdd.max <- g.p$sdd.max
  sdd.rate <- g.p$sdd.rate
  bird.hab <- g.p$bird.hab
  ncell <- nrow(sdd.alter)
  
  # initialize landscape & storage objects
  nbr <- 2 * sdd.max + 1
  sdd.new <- array(0, dim=c(nbr, nbr, ncell))
  bird.hab.E <- as.matrix(lc.df[,lc.col]) %*% (bird.hab %>% divide_by(sum(.)))
  
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
  
  for(n in 1:ncell) {
    ib <- sdd.0[,,2,sdd.alter$id.in[n]] != 0
    sdd.new[,,n][ib] <- d.pr[ib] * bird.hab.E[sdd.0[,,2,sdd.alter$id.in[n]][ib]]
  }
  sdd.new[,,] <- apply(sdd.new[,,], 3, function(x) x/sum(x))
  sdd.sparse.ls <- vector("list", ncell)
  for(n in 1:ncell) {
    sdd.sparse.ls[[n]] <- c(sdd.new[,,n])
    names(sdd.sparse.ls[[n]]) <- c(sdd.0[,,2,sdd.alter$id.in[n]])
    sdd.sparse.ls[[n]] <- sdd.sparse.ls[[n]][sdd.sparse.ls[[n]]>0]
  }
  sdd.sparse <- data.frame(i=rep(1:ncell, times=map_int(sdd.sparse.ls, length)),
                           j=unlist(map(sdd.sparse.ls, ~as.numeric(names(.)))),
                           pr=unlist(sdd.sparse.ls))
  sdd.sparse$j.idin <- lc.df$id.in[match(sdd.sparse$j, lc.df$id)]
  return(list(i=sdd.new, sp=sdd.sparse.ls, sp.df=sdd.sparse))
}




#'Short distance dispersal
#'
#'Calculate the number of viable seeds deposited within the short distance
#'dispersal neighborhood of each cell, accounting for the distance from the
#'source cell, bird habitat preference, the proportion of fruits eaten by birds,
#'and seed viability post-digestion.
#'@param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'  while \code{id.in} indexes only inbound cells
#'@param Fr Tibble of fruits produced in each cell output from
#'  \code{\link{make_fruits}}
#'@param gamma Scalar: mean number of seeds per fruit
#'@param p.c.E Vector of proportion of fruits eaten by birds from
#'  \code{\link{cell_E}}
#'@param s.c Proportion of viable seeds post-digestion
#'@param sdd.sp List of named vectors with SDD probabilities for each cell
#'  output from \code{\link{sdd_set_probs}}
#'@param sdd.rate Rate parameter for SDD exponential kernel
#'@param sdd.st \code{Logical} denoting whether to implement short distance
#'  dispersal stochastically
#'@param edges Character taking the value of one of: \code{"wall", "sink",
#'  "none"} where \code{"wall"} results in a dispersal probability of 0 for all
#'  out-of-bound cells with no populations modeled, \code{"sink"} results in
#'  dispersal of seeds to out-of-bound cells but no populations modeled, and
#'  \code{"none"} results in dispersal of seeds and populations modeled
#'@return Tibble with grid id and number of seeds in each cell
#'@keywords dispersal, SDD
#'@export

sdd_disperse <- function(id.i, Fr, gamma, p.c.E, s.c, 
                         sdd.sp, sdd.rate, sdd.st=T, edges="wall") {
  
  library(tidyverse); library(magrittr)
  
  # calculate seeds deposited within source cell vs emigrants
  N.source <- Fr %>%
    mutate(N.produced=round(gamma*N.fruit),
           N.emig=N.produced*p.c.E[id,]*s.c*(1-pexp(.5,sdd.rate)),
           N.dep=N.produced*p.c.E[id,]*s.c*pexp(.5,sdd.rate) + 
             N.produced*(1-p.c.E[id,])) %>%
    mutate(id.in=id.i$id.in[id]) %>%
    mutate_at(2:7, round)
  while(any(N.source$N.emig > .Machine$integer.max)) {
    ovr <- which(N.source$N.emig > .Machine$integer.max)
    for(i in seq_along(ovr)) {
      N.source$N.emig[ovr[i]] <- N.source$N.emig[ovr[i]] - .Machine$integer.max
      N.source <- rbind(N.source, (N.source[ovr[i],] %>%
                                       mutate(N.emig=.Machine$integer.max,
                                              N.dep=0, N.produced=0)))
    }
  }
  N.seed <- N.source %>% select(id, N.dep)
  
  if(sum(N.source$N.emig>0)>0) {
    if(sdd.st) {
      N.seed$N.dep <- round(N.seed$N.dep)
      SDD.sd <- vector("list", nrow(N.source))
      if(edges=="none") {
        for(n in 1:nrow(N.source)) {
          tmp <- rmultinom(1, N.source$N.emig[n], sdd.sp[[N.source$id[n]]])
          SDD.sd[[n]] <- tmp[tmp>0,1]
        }
      } else {
        for(n in 1:nrow(N.source)) {
          tmp <- rmultinom(1, N.source$N.emig[n], sdd.sp[[N.source$id.in[n]]])
          SDD.sd[[n]] <- tmp[tmp>0,1]
        }
      }
      SDD.dep <- unlist(SDD.sd)
      N.seed <- add_row(N.seed, 
                        id=as.numeric(names(SDD.dep)), 
                        N.dep=SDD.dep)
    } else {
      # assign emigrants to target cells & sum within each cell
      if(edges=="none") {
        N.seed <- N.seed %>% 
          add_row(id=unlist(apply(N.source, 1, function(x) names(sdd.sp[[x[1]]]))), 
                  N.dep=unlist(apply(N.source, 1, 
                                function(x) c(x[6] * sdd.sp[[x[1]]])))) %>%
          filter(N.dep > 0)
      } else {
        N.seed <- N.seed %>% 
          add_row(id=as.numeric(unlist(apply(N.source, 1, 
                                             function(x) names(sdd.sp[[x[8]]])))), 
                  N.dep=unlist(apply(N.source, 1, 
                                     function(x) c(x[6] * sdd.sp[[x[8]]])))) %>%
          filter(N.dep > 0)
      }
    }
  }
  
  N.seed <- N.seed %>% group_by(id) %>% 
    summarise(N=sum(N.dep) %>% round) %>%
    filter(N > 0)
  if(edges=="wall") N.seed <- filter(N.seed, !is.na(id.i$id.in[id]))
  
  return(list(N.seed=N.seed, N.source=N.source))
}




#' Long distance dispersal
#'
#' Assign n.ldd random long distance dispersal events across the landcape. A
#' single established seed is added to each target cell
#' @param ncell Number of inbound grid cells
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.in} indexes only inbound cells
#' @param M.0 \code{M.0} output from \code{\link{new_seedlings}} with grid id
#'   and number of new recruits in each cell
#' @param n.ldd Number of long distance dispersal events per time step
#' @return Tibble with grid id and number of seeds
#' @keywords LDD, dispersal
#' @export

ldd_disperse <- function(ncell, id.i, M.0, n.ldd) {
  
  if(n.ldd > 0) {
    ldd.id <- id.i$id[which(id.i$id.in %in% sample(1:ncell, n.ldd, replace=T))]
    M.0[ldd.id] <- M.0[ldd.id] + 1
  }

  return(M.0)
}







#' Short distance dispersal: lambda-based
#'
#' Calculate the number of immigrants to each cell, accounting for distance from
#' source cell and bird habitat preferences. It uses the output from simple
#' lambda-based population growth rather than the demographic version. The
#' output values account for carrying capacity.
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

sdd_lambda <- function(N.new, id.i, sdd.pr, sdd.rate, K.E, sdd.st=T) {
  # Calculate (N.arrivals | N.new, sdd.probs)
  # Accounts for distance from source cell & bird habitat preference
  # Returns dataframe with total population sizes.
  library(tidyverse)
  
  N.source <- N.new %>% filter(N.new > 0) %>% mutate(id.in=id.i$id.in[id])
  N.emig <- tibble(id=id.i$id, 
                   N=N.new$N.pop.upd[match(id.i$id, N.new$id)]) %>% 
    mutate(id.in=id.i$id.in[id])
  if(nrow(N.source) > 0) {
    if(sdd.st) {
      SDD.sd <- unlist(apply(N.source, 1,
                             function(x) sample(sdd.pr[,,2,x[5]], 
                                                round(x[3] * (1-pexp(.5,sdd.rate))), 
                                                replace=TRUE,
                                                prob=sdd.pr[,,1,x[5]])))
      SDD.dep <- tabulate(SDD.sd)  # vector of counts for 1:max(SDD.sd)
      SDD.nonzero <- SDD.dep > 0  # cell id's with N.dep > 0
      if(sum(SDD.nonzero)>0) {
        N.emig %<>% add_row(id=which(SDD.nonzero), 
                            N=SDD.dep[SDD.nonzero],
                            id.in=id.i$id.in[match(id, id.i$id)])
      }
    } else {
      N.emig %<>%
        add_row(id=apply(N.source, 1, function(x) c(sdd.pr[,,2,x[5]])) %>% c, 
                N=apply(N.source, 1, function(x) c(x[3] * (1-pexp(.5,sdd.rate)) * 
                                                     sdd.pr[,,1,x[5]])) %>% c,
                id.in=id.i$id.in[match(id, id.i$id)]) 
    }
  }
  N.emig %<>% filter(id != 0) %>% group_by(id) %>%
    summarise(N=sum(N, na.rm=T))
  N.emig$N <- round(pmin(K.E[N.emig$id,], N.emig$N))
  
  return(N.emig)
}




