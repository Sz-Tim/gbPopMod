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
