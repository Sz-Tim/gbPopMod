#' Assign landowner treatments
#'
#' This function randomly assigns pixels to the specified buckthorn management
#' treatment. Assigned cells either overwrite the cells treated in the previous
#' time step or are appended.
#' @param id_i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id_inbd} indexes only inbound cells
#' @param ncell Number of inbound grid cells
#' @param nTrt Proportion of inbound grid cells to be treated. If add != NULL,
#'   then this is appended to the cells treated in the previous time step.
#' @param trt.eff Named vector with effects of each treatment. These names are
#'   sampled and assigned to cells
#' @param addOwners \code{Logical} denoting whether to append new treated cells
#'   to the previous time step. If \code{TRUE}, then \code{trt.m1} must be
#'   provided
#' @param trt.m1 Tibble with grid id and treatment type from time \code{t - 1},
#'   either empty or generated with \code{\link{trt_assign}} in the previous
#'   time step.
#' @return Tibble with grid id and treatment type
#' @keywords control, treatment, owners
#' @export


trt_assign <- function(id_i, ncell, nTrt, trt.eff, 
                       addOwners=FALSE, trt.m1=NULL) {
  
  require(tidyverse)
  
  trt.t <- tibble(id=id_i$id[which(id_i$id_inbd %in% sample(1:ncell, nTrt))],
                  Trt=sample(names(trt.eff), nTrt, replace=TRUE))
  
  if(addOwners) {
    return(trt.m1 %<>% add_row(id=trt.t$id, Trt=trt.t$Trt))
  } else {
    return(trt.t)
  }
  
}




#' Implement ground cover treatments
#'
#' This function overwrites the establishment probabilities for cells receiving
#' a ground cover treatment
#' @param est.trt Tibble output from \code{\link{trt_assign}} with the grid id
#'   and treatment type for each cell
#' @param grd.trt Named vector with treatment types and associated establishment
#'   probability
#' @return Tibble with grid id and establishment probabilities
#' @keywords control, treatment, manual, litter, compaction, cover crop
#' @export

trt_ground <- function(est.trt, grd.trt) {
  
  require(tibble)
  
  pr.est.trt <- tibble(id=est.trt$id,
                       pr.est=grd.trt[match(est.trt$Trt, names(grd.trt))])
  
  return(pr.est.trt)
}




#' Implement cutting & spraying treatments
#'
#' This function adjusts the abundances within cells receiving cutting and/or
#' spraying treatments.
#' @param N.t Matrix or array of abundances, with dims=c(ngrid, (lc), y.ad)
#' @param y.ad Max(age at maturity)
#' @param N.trt Tibble output from \code{\link{trt_assign}} with the grid id and
#'   treatment type for each cell
#' @param man.trt Named vector with treatment types and associated success
#'   (=mortality) rates
#' @return Matrix or array of the same dimensions as N.t with adjusted
#'   abundances
#' @keywords control, treatment, manual, cutting, spraying
#' @export

trt_manual <- function(N.t, y.ad, N.trt, man.trt) {
  
  trt.eff <- tibble(id=N.trt$id,
                    surv=1-man.trt[match(N.trt$Trt, names(man.trt))])
  
  if(length(dim(N.t)) == 2) {
    N.t[trt.eff$id,] <- round(N.t[trt.eff$id,] * trt.eff$surv)
  } else {
    N.t[trt.eff$id,,] <- round(N.t[trt.eff$id,,] * trt.eff$surv)
  }
  
  return(N.t)
}




#' Change forested land cover to open invasible
#'
#' Specified pixels have a specified proportion of specified forest type 
#' converted to open invasible habitat instead. This function updates the land
#' cover proportions within \code{lc.df} and adjusts the short distance 
#' dispersal kernels in \code{sdd.pr} since bird habitat preferences will alter
#' the dispersal to each cell with land cover changes. This is a one way,
#' deterministic change with no natural reversion to forest. The affected 
#' parameters \emphasis{must} be updated with a call to \code{\link{cell_agg}}.
#' @param lc_chg.df Dataframe or tibble with cell grid indexes and a column for each forest type corresponding with the columns in \code{lc.df} specifying what proportion of forest should be converted to open invasible
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd.pr Array with \code{dim=c(i:disp.rows, j:disp.cols, k:2, n:ncell)} 
#'   output from \code{\link{sdd_set_probs}}
#' @param sdd.rate Rate parameter for SDD exponential kernel
#' @return Array with dim(disp.rows, disp.cols, 2, ncell) where the third
#'   dimension contains grid id's for the neighborhood or probabilities to each
#'   target cell
#' @keywords sdd, dispersal, probability, probabilities
#' @export

clear_forest <- function(lc_ch.df, lc.df, sdd.pr, sdd.rate) {
  # this happens
  return(list(lc.df=lc.df, sdd.pr=sdd.pr))
}
