#' Initialize populations randomly
#'
#' This function initializes populations randomly with populated cells
#' containing adults at 50% K and juveniles at 10% K
#' @param ngrid Number of grid cells in entire map
#' @param g.p Named list of global parameters
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @return Matrix or array of initial abundances with dim=c(ngrid, (n.lc), y.ad)
#' @keywords initialize, set up
#' @export


pop_init <- function(ngrid, g.p, lc.df) {
  
  p.0 <- sample(lc.df$id[lc.df$inbd], g.p$N.p.t0)
  y.ad <- max(g.p$age.f)  # adult age bin
  
  if(length(g.p$age.f) == 1) {
    N.init <- matrix(0, ngrid, y.ad)  # column for each age class
    N.init[p.0,y.ad] <- round(as.matrix(lc.df[lc.df$id %in% p.0,4:9]) %*% 
                                (g.p$K/2))
    N.init[p.0,-y.ad] <- round(N.init[p.0,y.ad]/5)
    
  } else {
    N.init <- array(0, dim=c(ncell, g.p$n.lc, y.ad))
    N.init[p.0,,y.ad] <- round(t(t(as.matrix(lc.df[lc.df$id %in% p.0,4:9])) * 
                                   g.p$K/2))
    N.init[p.0,,-y.ad] <- round(N.init[p.0,,y.ad]/5)
  }
  
  return(N.init)
}