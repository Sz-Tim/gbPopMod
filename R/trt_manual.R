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
  # given cells, treatments, and treatment effects: adjust N.adults
  
  # spraying probably has a %kill for all N
  # some % die, some % get bumped back to juvenile
  # effects on germination rates?
  
  prop.kill <- tibble(id=N.trt$id,
                      prop=1-man.trt[match(N.trt$Trt, names(man.trt))])
  if(length(dim(N.t)) == 2) {
    N.t[prop.kill$id,] <- round(N.t[prop.kill$id,] * prop.kill$prop)
  } else {
    N.t[prop.kill$id,,] <- round(N.t[prop.kill$id,,] * prop.kill$prop)
  }
  return(N.t)
}