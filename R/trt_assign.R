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