#' Assign landowner treatments
#'
#' This function randomly assigns pixels to the specified buckthorn management
#' treatment. Assigned cells either overwrite the cells treated in the previous
#' time step or are appended.
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.inbd} indexes only inbound cells
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


trt_assign <- function(id.i, ncell, nTrt, trt.eff, 
                       addOwners=FALSE, trt.m1=NULL) {
  
  library(tidyverse)
  
  trt.t <- tibble(id=id.i$id[which(id.i$id.inbd %in% sample(1:ncell, nTrt))],
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
  
  library(tibble)
  
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




#' Assign cells to convert a proportion of forest to open invasible
#'
#' Randomly assigns a specified number of cells to have their forest habitat
#' cleared by a random proportion.
#' @param nChg Proportion of inbound cells to cut
#' @param ncell Number of inbound grid cells
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param f.c Vector of forest column indexes within \code{lc.df}
#' @return Tibble \code{id.chg} with grid and inbound indexes of cells with land
#'   cover change and matrix \code{mx} with the change in each forest category
#'   and the total change to be added to open habitat
#' @keywords control, cut, forest, land cover change, owners
#' @export

cut_assign <- function(nChg, ncell, lc.df, f.c) {
  
  library(tidyverse)
  
  id.lc <- sample(1:ncell, nChg*ncell)
  n <- length(id.lc)
  id.chg <- dplyr::filter(lc.df, id.inbd %in% id.lc) %>% 
    select(id, id.inbd)
  mx <- (runif(n*length(f.c)) * lc.df[id.chg$id, f.c]) %>%
    cbind(., TotChg=rowSums(.))
  
  return(list(id.chg=id.chg, mx=mx))
}





#' Change forested land cover to open invasible
#'
#' Specified pixels have a specified proportion of specified forest type
#' converted to open invasible habitat. This function updates the land cover
#' proportions within \code{lc.df}. This is a one way, deterministic change with
#' no natural reversion to forest. The affected parameters \emph{must} be
#' updated with a call to \code{\link{cell_agg}} and to
#' \code{\link{sdd_set_probs}(lc.new=id.chg)}.
#' @param id.chg Dataframe or tibble with cell grid indexes identifying which
#'   cells are to change. if not input manually, can be randomly created by
#'   \code{\link{cut_assign}}
#' @param f.chg Matrix with a column for each forest type corresponding with the
#'   forest columns in \code{lc.df} specifying the amount of each category of
#'   forest to be converted, and a final column \code{TotChg} with the total
#'   amount converted to open invasible. If not input manually, can be randomly
#'   created by \code{\link{cut_assign}}
#' @param f.c Vector of forest column indexes within \code{lc.df}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @return List with \code{lc.df} and \code{sdd.pr} where both are sparse,
#'   containing only updated values for the cells with altered land cover.
#' @keywords control, cut, forest, land cover change, owners
#' @export

cut_forest <- function(id.chg, f.chg, f.c, lc.df) {
  
  library(tidyverse)
  
  # shift forest to open 
  chg.df <- lc.df[id.chg$id,]
  chg.df[,f.c] <- chg.df[,f.c] - f.chg[,1:length(f.c)]
  chg.df[,"Opn"] <- chg.df[,"Opn"] + f.chg[,length(f.c)+1]

  return(chg.df)
}
