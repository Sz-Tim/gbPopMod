#' Assign landowner treatments
#'
#' Assign pixels to the specified buckthorn management treatment. Assigned cells
#' either overwrite the cells treated in the previous time step or are appended.
#' @param id.i Tibble matching cell IDs. \code{id} indexes on the entire grid
#'   while \code{id.in} indexes only inbound cells
#' @param ncell \code{NULL} Number of inbound grid cells. Required if cell IDs
#'   are not supplied via \code{assign_i}
#' @param assign_i \code{NULL} Vector of inbound cell IDs to treat. If
#'   \code{NULL}, then \code{ncell} IDs are sampled randomly for treatments
#' @param pTrt Proportion of inbound grid cells to be treated. If \code{add !=
#'   NULL}, then this is appended to the cells treated in the previous time
#'   step.
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

trt_assign <- function(id.i, ncell=NULL, assign_i=NULL, pTrt, trt.eff, 
                       addOwners=FALSE, trt.m1=NULL) {
  
  library(tidyverse)
  
  if(is.null(assign_i)) {
    assign_i <- sample(1:ncell, ceiling(pTrt*ncell))
  } else {
    nTrt <- length(assign_i)
  }
  
  trt.t <- tibble(id=id.i$id[which(id.i$id.in %in% assign_i)],
                  Trt=sample(names(trt.eff), ceiling(pTrt*ncell), replace=TRUE))
  
  if(addOwners) {
    return(trt.m1 %<>% add_row(id=trt.t$id, Trt=trt.t$Trt))
  } else {
    return(trt.t)
  }
  
}




#' Implement ground cover treatments
#'
#' Overwrite the establishment probabilities for cells receiving a ground cover
#' treatment
#' @param est.trt Tibble output from \code{\link{trt_assign}} with the grid id
#'   and treatment type for each cell
#' @param grd.trt Named vector with treatment types and associated establishment
#'   probability
#' @return Tibble with grid id and establishment probabilities
#' @keywords control, treatment, manual, litter, compaction, cover crop
#' @export

trt_ground <- function(est.trt, grd.trt) {
  
  library(tibble)
  
  p.trt <- tibble(id=est.trt$id,
                       p=grd.trt[match(est.trt$Trt, names(grd.trt))])
  
  return(p.trt)
}




#' Implement cutting & spraying treatments
#'
#' Adjust the abundances within cells receiving cutting and/or spraying
#' treatments.
#' @param N.t Matrix or array of abundances, with dims=c(ngrid, (lc), m.max)
#' @param m.max Max(age at maturity)
#' @param N.trt Tibble output from \code{\link{trt_assign}} with the grid id and
#'   treatment type for each cell
#' @param man.trt Named vector with treatment types and associated success
#'   (=mortality) rates
#' @return Matrix or array of the same dimensions as N.t with adjusted
#'   abundances
#' @keywords control, treatment, manual, cutting, spraying
#' @export

trt_manual <- function(N.t, m.max, N.trt, man.trt) {
  
  library(tibble)
  
  trt.eff <- tibble(id=N.trt$id,
                    surv=1-man.trt[match(N.trt$Trt, names(man.trt))])
  
  if(length(dim(N.t)) == 2) {
    N.t[trt.eff$id,] <- round(N.t[trt.eff$id,] * trt.eff$surv)
  } else {
    N.t[trt.eff$id,,] <- round(N.t[trt.eff$id,,] * trt.eff$surv)
  }
  
  return(N.t)
}




#' Assign cells to convert a proportion of forest.col to open invasible
#'
#' Assign a specified number of random cells to have their forest.col habitat
#' cleared by a random proportion.
#' @param pChg Proportion of inbound cells to cut
#' @param ncell \code{NULL} Number of inbound grid cells. Required if cell IDs
#'   are not supplied via \code{assign_i}
#' @param assign_i \code{NULL} Vector of inbound cell IDs to treat. If
#'   \code{NULL}, then \code{ncell} IDs are sampled randomly for treatments
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param forest.col Vector of forest.col column indexes within \code{lc.df}
#' @return Tibble \code{id.chg} with grid and inbound indexes of cells with land
#'   cover change and matrix \code{mx} with the change in each forest.col
#'   category and the total change to be added to open habitat
#' @keywords control, cut, forest.col, land cover change, owners
#' @export

cut_assign <- function(pChg, ncell=NULL, assign_i=NULL, lc.df, forest.col) {
  
  library(tidyverse)
  
  if(is.null(assign_i)) assign_i <- sample(1:ncell, ceiling(pChg*ncell))
  n <- length(assign_i)
  id.chg <- dplyr::filter(lc.df, id.in %in% assign_i) %>% 
    select(id, id.in)
  mx <- (runif(n*length(forest.col)) * lc.df[id.chg$id, forest.col]) %>%
    cbind(., TotChg=rowSums(.))
  
  return(list(id.chg=id.chg, mx=mx))
}





#' Change forested land cover to open invasible
#'
#' Convert forest to open habitat. The specified pixels have a specified
#' proportion of specified forest.col type converted to open invasible habitat.
#' This function updates the land cover proportions within \code{lc.df}. This is
#' a one way, deterministic change with no natural reversion to forest.col. The
#' affected parameters \emph{must} be updated with a call to
#' \code{\link{cell_agg}} and to \code{\link{sdd_set_probs}(lc.new=id.chg)}.
#' @param id.chg Dataframe or tibble with cell grid indexes identifying which
#'   cells are to change. if not input manually, can be randomly created by
#'   \code{\link{cut_assign}}
#' @param forest.chg Matrix with a column for each forest.col type corresponding
#'   with the forest.col columns in \code{lc.df} specifying the amount of each
#'   category of forest.col to be converted, and a final column \code{TotChg}
#'   with the total amount converted to open invasible. If not input manually,
#'   can be randomly created by \code{\link{cut_assign}}
#' @param forest.col Vector of forest.col column indexes within \code{lc.df}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @return List with \code{lc.df} and \code{sdd.pr} where both are sparse,
#'   containing only updated values for the cells with altered land cover.
#' @keywords control, cut, forest.col, land cover change, owners
#' @export

cut_forest <- function(id.chg, forest.chg, forest.col, lc.df) {
  
  library(tidyverse)
  
  # shift forest.col to open 
  chg.df <- lc.df[id.chg$id,]
  chg.df[,forest.col] <- chg.df[,forest.col] - forest.chg[,1:length(forest.col)]
  chg.df[,"Opn"] <- chg.df[,"Opn"] + forest.chg[,length(forest.col)+1]

  return(chg.df)
}
