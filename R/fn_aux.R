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




#' Expand to all combinations of land cover parameter ranges
#'
#' This function is similar to \code{\link[base]{expand.grid}} but inputs six
#' vectors with min and max parameter values (one set for each land cover
#' category) in addition to a length.out parameter, and returns a list of
#' vectors, with an element for each combination of land cover parameters.
#' @param lc.min Vector \code{length=n.lc} of minimum parameter values
#' @param lc.max Vector \code{length=n.lc} of maximum parameter values
#' @param len.out Number of parameter values per land cover category
#' @return List of vectors, each length 6
#' @keywords expand.grid, sensitivity
#' @export

expand_all_LCs <- function(lc.min=rep(0.1, 6), lc.max=rep(0.9, 6), len_out=2) {
  library(tidyverse); library(purrr)
  names(lc.min) <- c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")
  names(lc.max) <- names(lc.min)
  g <- map2(lc.min, lc.max, ~seq(.x, .y, length.out=len_out)) %>% 
    expand.grid %>% as.matrix
  return(l=lapply(seq_len(nrow(g)), function(i) g[i,]))
}




#' Expand to all combinations of canopy parameter ranges
#'
#' This function is similar to \code{\link[base]{expand.grid}} but inputs two
#' vectors with min and max parameter values (one set for open canopy and one
#' set for closed canopy) in addition to a length.out parameter, and returns a
#' list of vectors, with an element for each combination of land cover
#' parameters.
#' @param Open Vector of min & max for open canopy categories (open invasible
#'   and other)
#' @param Closed Vector of min & max for closed canopy categories (deciduous,
#'   white pine, other evergreen, and mixed forests)
#' @param length.out Number of parameter values per land cover category
#' @return List of vectors, each length 6
#' @keywords expand.grid, sensitivity
#' @export

expand_cnpy <- function(Op=c(3, 7), Cl=c(3, 7), length_out=2) {
  g <- expand.grid(Op=seq(Op[1], Op[2], length.out=length_out),
                   Cl=seq(Cl[1], Cl[2], length.out=length_out))
  g <- t(apply(g, 1, function(i) rep(i, times=c(2,4))))
  return(l=lapply(seq_len(nrow(g)), function(i) g[i,]))
}




#' Aggregate compositional data within each cell
#'
#' This function reformats and calculates cell-means based on land cover
#' composition for relevant parameters.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param K Vector \code{length=n.lc} with carrying capacity for each land cover
#'   type
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
#' @param edges Character taking the value of one of: \code{"wall", "sink",
#'   "none"} where \code{"wall"} results in a dispersal probability of 0 for all
#'   out-of-bound cells with no populations modeled, \code{"sink"} results in
#'   dispersal of seeds to out-of-bound cells but no populations modeled, and
#'   \code{"none"} results in dispersal of seeds and populations modeled
#' @return Named list with values aggregated within cells based on land cover
#'   types. Includes: \describe{ \item{\code{lc.mx}}{Matrix \code{(ncol=n.lc,
#'   nrow=ngrid)} with land cover proportions} \item{\code{K.ag}}{Vector
#'   \code{length=ngrid} with total K} \item{\code{K.lc}}{Matrix
#'   \code{(ncol=n.lc, nrow=ngrid)} with K per land cover category}
#'   \item{\code{pr.s.ag}}{Vector \code{length=ngrid} with pr(surv)}
#'   \item{\code{rel.dens}}{Matrix \code{(ncol=n.lc, nrow=ngrid)} with relative
#'   density among land cover categories} \item{\code{fec.ag}}{Vector
#'   \code{length=ngrid} with mean fruit produced per adult)}
#'   \item{\code{pr.f.ag}}{Vector \code{length=ngrid} with fruiting probability}
#'   \item{\code{pr.eat.ag}}{Vector \code{length=ngrid} with proportion eaten by
#'   birds} \item{\code{pr.est.ag}}{Vector \code{length=ngrid} with seedling
#'   establishment probabilities} }
#' @note If \code{!is.null(pr.est.trt)}, then the associated pr.est.ag values
#'   are substituted in the cells that received a relevant management
#'   treatments.
#' @keywords premultiply, aggregate, set up, initialize
#' @export

cell_agg <- function(lc.df, K, pr.s, fec, pr.f, pr.eat, 
                     pr.est, pr.est.trt=NULL, edges="wall") {
  
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
  
  if(edges=="sink") pr.est.ag[!lc.df$inbd] <- 0
  
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
#' @return Matrix or array of initial abundances with \code{dim=c(ngrid, (n.lc),
#'   y.ad)}
#' @keywords initialize, set up
#' @export

pop_init <- function(ngrid, g.p, lc.df) {
  
  p.0 <- sample(lc.df$id[lc.df$inbd], g.p$N.p.t0)
  y.ad <- max(g.p$age.f)  # adult age bin
  
  if(length(g.p$age.f) == 1) {
    N.init <- matrix(0, ngrid, y.ad)  # column for each age class
    N.init[p.0,y.ad] <- round(as.matrix(lc.df[lc.df$id %in% p.0, 4:9]) %*% 
                                (g.p$K/2))
    N.init[p.0,-y.ad] <- round(N.init[p.0,y.ad]/5)
    
  } else {
    N.init <- array(0, dim=c(ngrid, g.p$n.lc, y.ad))
    N.init[p.0,,y.ad] <- round(t(t(as.matrix(lc.df[lc.df$id %in% p.0, 4:9])) * 
                                   g.p$K/2))
    N.init[p.0,,-y.ad] <- round(N.init[p.0,,y.ad]/5)
  }
  
  return(N.init)
}




#' Set global parameters
#'
#' This function sets the global parameters for the simulation, allowing for
#' individual elements to be reassigned.
#' @param tmax \code{100} Number of time steps per simulation
#' @param dem.st \code{FALSE} Include stochasticity in demography?
#' @param sdd.st \code{TRUE} Include stochasticity in short distance dispersal?
#' @param bank \code{TRUE} Include seedbank?
#' @param n.cores \code{4} Number of cores for parallelizing sdd.pr calculation
#' @param lc.r \code{100} Maximum number of rows (\code{y}) in landscape
#' @param lc.c \code{100} Maximum number of columns (\code{x}) in landscape
#' @param n.lc \code{6} Number of land cover categories
#' @param N.p.t0 \code{10} Number of cells with buckthorn at t=1
#' @param K \code{c(750, 10, 100, 100, 300, 100)} Vector (length=n.lc) of
#'   carrying capacities for adults
#' @param pr.s \code{c(0.9, 0.1. 0.6, 0.6, 0.6, 0.6)} Vector \code{length=n.lc}
#'   of annual juvenile survival rates
#' @param pr.f \code{c(0.9, 0.1, 0.29, 0.23, 0.2, 0.3)} Vector
#'   \code{length=n.lc} of fruiting probabilities
#' @param fec \code{c(200, 100, 40, 20, 20, 10)} Vector \code{length=n.lc} of
#'   mean fruit per adult
#' @param age.f \code{4} Vector \code{length=n.lc} or scalar of age at first
#'   fruiting. Individuals at this age are considered adults
#' @param pr.sb \code{0.3} Probability of annual survival in seed bank
#' @param pr.est \code{c(0.07, 0.01, 0.08, 0.02, 0.02, 0.03)} Vector
#'   \code{length=n.lc} of seedling establishment probabilities
#' @param sdd.max \code{15} Maximum dispersal distance in cells
#' @param sdd.rate \code{0.1} 1/mn for exponential dispersal kernel
#' @param n.ldd \code{1} Number of long distance dispersal events per year
#' @param pr.eat \code{c(0.3, 0.1, 0.2, 0.2, 0.2, 0.1)} Vector
#'   \code{length=n.lc} of proportion of fruits eaten by birds, with
#'   \code{1-pr.eat} assumed to drop directly below buckthorn individuals
#' @param bird.hab \code{c(0.35, 0.35, 0.05, 0.1, 0.1, 0.05)} Vector
#'   \code{length=n.lc} of bird habitat preferences
#' @param pr.s.bird \code{0.6} Seed viability post-digestion
#' @param edges \code{wall} Boundary behavior, taking values of \code{wall},
#'   \code{sink}, or \code{none}. See boundary_behavior.Rmd for descriptions
#' @return Named list of global parameters including all arguments as elements
#' @keywords initialize, set up, global, parameter
#' @export

set_g_p <- function(tmax=100, dem.st=FALSE, sdd.st=TRUE, bank=TRUE, n.cores=4, 
                    lc.r=100, lc.c=100, n.lc=6, N.p.t0=10,
                    K=c(750, 10, 100, 100, 300, 100),
                    pr.s=c(0.9, 0.1, 0.6, 0.6, 0.6, 0.6),
                    pr.f=c(0.9, 0.1, 0.29, 0.23, 0.2, 0.3),
                    fec=c(200, 100, 40, 20, 20, 10),
                    age.f=4, pr.sb=0.3, 
                    pr.est=c(0.07, 0.01, 0.08, 0.02, 0.02, 0.03),
                    sdd.max=15, sdd.rate=0.1, n.ldd=1,
                    pr.eat=c(0.3, 0.1, 0.2, 0.2, 0.2, 0.1),
                    bird.hab=c(.35, .35, 0.05, 0.1, 0.1, 0.05), pr.s.bird=0.6,
                    edges="wall") {
  
  g.p <- list(tmax=tmax, dem.st=dem.st, sdd.st=sdd.st, bank=bank,
              n.cores=n.cores, lc.r=lc.r, lc.c=lc.c, n.lc=n.lc, N.p.t0=N.p.t0, 
              K=K, pr.s=pr.s, pr.f=pr.f, fec=fec, age.f=age.f, pr.sb=pr.sb,
              pr.est=pr.est, sdd.max=sdd.max, sdd.rate=sdd.rate, n.ldd=n.ldd,
              pr.eat=pr.eat, bird.hab=bird.hab, pr.s.bird=pr.s.bird, 
              edges=edges)
  
  return(g.p)
}




#' Set buckthorn control treatment parameters
#'
#' This function sets the buckthorn control treatment parameters for the
#' simulation, allowing for individual elements to be reassigned.
#' @param null_ctrl \code{TRUE} Set control parameters to \code{NULL}?
#' @param t.trt \code{30} Year to start treatments
#' @param add.owners \code{FALSE} Do owners treat every year once starting a
#'   particular treatment?
#' @param nTrt.grd \code{0.05} Proportion of cells with ground treatments in
#'   each time step
#' @param nTrt.man \code{0.05} Proportion of cells with manual treatments in
#'   each time step
#' @param grd.trt \code{Lit=0.005, Cov=0.01, Com=0.00001} Named vector with
#'   ground treatments and associated seedling establishment probabilities
#' @param man.trt \code{c(M=0.1, C=0.3, MC=0.8)} Named vector with manual
#'   treatments and associated mortality (=success) rates
#' @param lc.chg \code{TRUE} Does land cover change across years?
#' @param n.chg \code{0.0001} Proportion of cells with land cover change each
#'   year
#' @return Named list of control parameters including all arguments as elements
#'   unless \code{null_ctrl==TRUE}, in which case the function returns
#'   \code{NULL}
#' @keywords initialize, set up, control, treatment, parameter
#' @export

set_control_p <- function(null_ctrl=TRUE, t.trt=30, add.owners=FALSE,
                          nTrt.grd=0.05, nTrt.man=0.05,
                          grd.trt=c(Lit=0.005, Cov=0.01, Com=0.00001),
                          man.trt=c(M=0.1, C=0.3, MC=0.8),
                          lc.chg=TRUE, n.chg=0.0001) {
  
  if(null_ctrl) {
    control.p <- NULL
  } else {
    control.p <- list(t.trt=t.trt, add.owners=add.owners,
                      nTrt.grd=nTrt.grd, nTrt.man=nTrt.man,
                      grd.trt=grd.trt, man.trt=man.trt,
                      lc.chg=lc.chg, n.chg=n.chg)
  }
  
  return(control.p)
}
