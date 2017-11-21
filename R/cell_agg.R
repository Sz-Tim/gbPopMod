#' Aggregate compositional data within each cell
#'
#' This function reformats and calculates cell-means based on land cover
#' composition for relevant parameters.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param K Vector with carrying capacity for each land cover type
#' @param pr.s Vector with juvenile survival probability for each land cover
#'   type
#' @param fec Vector with mean per-individual fruit production for each land
#'   cover type
#' @param pr.f Vector with mean probability of fruiting for each land cover type
#' @param pr.eat Vector with proportion of fruits eaten by birds for each land
#'   cover type
#' @param pr.est Vector with probability of establishment for each land cover
#'   type
#' @param pr.est.trt Tibble with grid id and modified establishment
#'   probabilities for cells with ground cover treatments; default = NULL
#' @return Named list with values aggregated within cells based on land cover
#'   types. Includes: \describe{ \item{lc.mx}{matrix(col=n.lc, row=ngrid) with
#'   LC proportions} \item{K.ag}{vector(length=ngrid) with total K}
#'   \item{K.lc}{matrix(col=n.lc, row=ngrid) with K per LC}
#'   \item{pr.s.ag}{vector(length=ngrid) with pr(surv)}
#'   \item{rel.dens}{matrix(col=n.lc, row=ngrid) with relative density among LC}
#'   \item{fec.ag}{vector(length=ngrid) with mn(fruit per adult)}
#'   \item{pr.f.ag}{vector(length=ngrid) with pr(fruit)}
#'   \item{pr.eat.ag}{vector(length=ngrid) with pr(eaten by bird)}
#'   \item{pr.est.ag}{vector(length=ngrid) with pr(establish)} }
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
