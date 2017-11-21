#' Long distance dispersal
#'
#' This function assigns n.ldd random long distance dispersal events across the
#' landcape. A single seed is added to each target cell
#' @param ncell Number of inbound grid cells
#' @param id_i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id_inbd} indexes only inbound cells
#' @param N.seed \code{N.seed} output from \code{\link{sdd_disperse}} with grid
#'   id and number of seeds in each cell
#' @param n.ldd Number of long distance dispersal events per time step
#' @return Tibble with grid id and number of seeds
#' @keywords LDD, dispersal
#' @export


ldd_disperse <- function(ncell, id_i, N.seed, n.ldd) {

  require(tidyverse); require(magrittr)
    
  ldd.id <- id_i$id[which(id_i$id_inbd == sample(1:ncell, n.ldd, replace=T))]

  N.seed %<>% 
    add_row(id=ldd.id, N=rep(1, n.ldd)) %>%
    group_by(id) %>% summarise(N=sum(N))
  
  return(N.seed)
}