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