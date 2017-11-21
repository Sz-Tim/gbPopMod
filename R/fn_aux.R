#' Expand all pairwise combinations of two vectors into one character vector
#'
#' This function is similar to \code{\link[base]{expand.grid}} but inputs two
#' vectors and returns a single character vector with the values from the two
#' vectors separated by "_" by default.
#' @param x Vector
#' @param y Vector
#' @param sep Like paste, specifies character(s) to separate vector values
#' @return Vector (character)
#' @keywords expand.grid
#' @export


expand_v <- function(x, y, sep="_") {
  paste(rep.int(x, length(y)), 
        rep.int(y, rep.int(length(x),length(y))),
        sep=sep)
}




#' Aggregate compositional data within each cell
#'
#' This function reformats and calculates cell-means based on land cover
#' composition for relevant parameters.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param K Vector \code{length=n.lc} with carrying capacity for each land 
#'   cover type
#' @param pr.s Vector \code{length=n.lc} with juvenile survival probability for
#'   each land cover type
#' @param fec Vector \code{length=n.lc} with mean per-individual fruit 
#'   production for each land cover type
#' @param pr.f Vector \code{length=n.lc} with mean probability of fruiting for 
#'   each land cover type
#' @param pr.eat Vector \code{length=n.lc} with proportion of fruits eaten by 
#'   birds for each land cover type
#' @param pr.est Vector \code{length=n.lc} with seedling establishment 
#'   probability for each land cover type
#' @param pr.est.trt Tibble with grid id and modified establishment
#'   probabilities for cells with ground cover treatments; default = NULL
#' @return Named list with values aggregated within cells based on land cover
#'   types. Includes: 
#'   \describe{ 
#'     \item{\code{lc.mx}}{Matrix \code{(ncol=n.lc, nrow=ngrid)} with land
#'       cover proportions} 
#'     \item{\code{K.ag}}{Vector \code{length=ngrid} with total K}
#'     \item{\code{K.lc}}{Matrix \code{(ncol=n.lc, nrow=ngrid)} with K per land
#'       cover category}
#'     \item{\code{pr.s.ag}}{Vector \code{length=ngrid} with pr(surv)}
#'     \item{\code{rel.dens}}{Matrix \code{(ncol=n.lc, nrow=ngrid)} with 
#'       relative density among land cover categories}
#'     \item{\code{fec.ag}}{Vector \code{length=ngrid} with mean fruit produced
#'       per adult)}
#'     \item{\code{pr.f.ag}}{Vector \code{length=ngrid} with fruiting 
#'       probability}
#'     \item{\code{pr.eat.ag}}{Vector \code{length=ngrid} with proportion
#'       eaten by birds}
#'     \item{\code{pr.est.ag}}{Vector \code{length=ngrid} with seedling
#'       establishment probabilities} 
#'   }
#' @note If \code{!is.null(pr.est.trt)}, then the associated pr.est.ag values
#'   are substituted in the cells that received a relevant management
#'   treatments.
#' @keywords premultiply, aggregate, set up, initialize
#' @export

cell_agg <- function(lc.df, K, pr.s, fec, pr.f, pr.eat, 
                     pr.est, pr.est.trt=NULL) {
  
  lc.mx <- as.matrix(lc.df[,4:9])
  K.ag <- round(lc.mx %*% K)
  K.lc <- round(t(t(lc.mx) * K))
  rel.dens <- t(apply(lc.mx, 1, function(x) K*x/c(x%*%K)))
  pr.s.ag <- c(lc.mx %*% pr.s)
  fec.ag <- lc.mx %*% fec
  pr.f.ag <- lc.mx %*% pr.f
  pr.eat.ag <- lc.mx %*% pr.eat
  pr.est.ag <- lc.mx %*% pr.est
  
  if(!is.null(pr.est.trt)) {
    pr.est.ag[pr.est.trt$id,] <- pr.est.trt$pr.est
  }
  
  return(list(lc.mx=lc.mx, K.ag=K.ag, K.lc=K.lc, rel.dens=rel.dens,
              pr.s.ag=pr.s.ag, fec.ag=fec.ag, pr.f.ag=pr.f.ag,
              pr.eat.ag=pr.eat.ag, pr.est.ag=pr.est.ag))
}




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