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