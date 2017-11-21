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
#' @param control.p \code{NULL} or named list of buckthorn control treatment
#'   parameters
#' @return Array with dim(disp.rows, disp.cols, 2, ncell) where the third
#'   dimension contains grid id's for the neighborhood or probabilities to each
#'   target cell
#' @keywords sdd, dispersal, probability, probabilities
#' @export

sdd_set_probs <- function(ncell, lc.df, g.p) {

  require(purrr); require(tidyverse); require(pbapply); require(fastmatch)
  
  # unpack parameters
  sdd.max <- g.p$sdd.max
  sdd.rate <- g.p$sdd.rate
  bird.hab <- g.p$bird.hab
  
  # initialize landscape & storage objects
  n.x <- range(lc.df[,1])
  n.y <- range(lc.df[,2])
  nbr <- 2 * sdd.max + 1
  sdd.i <- array(0, dim=c(nbr, nbr, 2, ncell))
  bird.hab.ag <- as.matrix(lc.df[,4:9]) %*% (bird.hab %>% divide_by(sum(.)))
  
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
  xx <- map(lc.df$x[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
  yy <- map(lc.df$y[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
  
  # create lists of on-map xy neighborhood ranges
  n_ix <- map(xx, ~.[.>=n.x[1] & .<=n.x[2]])
  n_iy <- map(yy, ~.[.>=n.y[1] & .<=n.y[2]])
  
  # generate all xy combinations & neighborhood matrix indices
  cat("generating neighborhoods...\n")
  n_i <- map2(n_ix, n_iy, expand_v) 
  n_x <- map2(xx, n_ix, `%fin%`) %>% map(which) %>% map(range)
  n_y <- map2(yy, n_iy, `%fin%`) %>% map(which) %>% map(range)
  
  # match xy combinations with cell IDs
  cat("determining neighborhood cell IDs...\n")
  pboptions(type="none")
  if(g.p$n_cores > 1) {
    p.c <- makeCluster(g.p$n_cores)
    c_i <- pblapply(n_i, function(x) fastmatch::fmatch(x, lc.df$x_y), cl=p.c)
    stopCluster(p.c)
  } else {
    c_i <- pblapply(n_i, function(x) fastmatch::fmatch(x, lc.df$x_y))
  }
  
  cat("calculating probabilities...\n")
  for(n in 1:ncell) {
    # find cell ID for each cell in neighborhood
    sdd.i[n_y[[n]][1]:n_y[[n]][2],
          n_x[[n]][1]:n_x[[n]][2],2,n] <- matrix(c_i[[n]], 
                                                 ncol=diff(n_x[[n]])+1,
                                                 byrow=TRUE)
    # weight by bird habitat preference
    ib <- sdd.i[,,2,n] != 0  # inbounds neighbors
    sdd.i[,,1,n][ib] <- d.pr[ib] * bird.hab.ag[sdd.i[,,2,n][ib]]
    # set cell ID to 0 if pr(target) == 0 
    sdd.i[,,2,n][sdd.i[,,1,n]==0] <- 0
    # progress update
    if(n %% 5000 == 0) {
      cat("finished cell", n, "\n")
    }
  }
  cat("finished:", n, "cells\n")
  sdd.i[,,1,] <- apply(sdd.i[,,1,], 3, function(x) x/sum(x))
  return(sdd.i)
}
