#' Run simulation
#'
#' This function runs the simulation. Currently, it runs all time steps, but for
#' the economic model, the structure will need to be slightly adjusted to run a
#' single time step. The initialization is separated from this function for that
#' reason.
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters including:
#'   \describe{
#'     \item{\code{tmax}}{Number of time steps per simulation}
#'     \item{\code{dem.st}}{\code{Logical} Include stochasticity in demography?}
#'     \item{\code{sdd.st}}{\code{Logical} Include stochasticity in short
#'       distance dispersal?}
#'     \item{\code{n_cores}}{Number of cores for parallelizing sdd.pr 
#'       calculation}
#'     \item{\code{lc.r}}{Maximum number of rows (\code{y}) in landscape}
#'     \item{\code{lc.c}}{Maximum number of columns (\code{x}) in landscape}
#'     \item{\code{n.lc}}{Number of land cover categories} 
#'     \item{\code{N.p.t0}}{Number of cells with buckthorn at t=1} 
#'     \item{\code{K}}{Vector (length=n.lc) of carrying capacities for adults} 
#'     \item{\code{pr.s}}{Vector \code{length=n.lc} of annual juvenile survival
#'       rates}
#'     \item{\code{pr.f}}{Vector \code{length=n.lc} of fruiting probabilities} 
#'     \item{\code{fec}}{Vector \code{length=n.lc} of mean fruit per adult} 
#'     \item{\code{age.f}}{Vector \code{length=n.lc} or scalar of age at first
#'       fruiting. Individuals at this age are considered adults}
#'      \item{\code{bank}}{\code{Logical} Include seedbank?} 
#'     \item{\code{pr.sb}}{Probability of annual survival in seed bank} 
#'     \item{\code{pr.est}}{Vector \code{length=n.lc} of seedling establishment 
#'       probabilities }
#'     \item{\code{sdd.max}}{Maximum dispersal distance in cells} 
#'     \item{\code{sdd.rate}}{1/mn for exponential dispersal kernel}
#'     \item{\code{pr.eat}}{Vector \code{length=n.lc} of proportion of fruits 
#'       eaten by birds, with \code{1-pr.eat} assumed to drop directly below
#'       buckthorn individuals}
#'     \item{\code{bird.hab}}{Vector \code{length=n.lc} of bird habitat 
#'       preferences}
#'     \item{\code{pr.s.bird}}{Seed viability post-digestion} 
#'     \item{\code{n.ldd}}{Number of long distance dispersal events per year} 
#'   }
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd.pr Array with sdd probabilities and neighborhoods created by 
#'   \code{\link{sdd_set_probs}}
#' @param N.init Matrix or array with initial population sizes created by 
#'   \code{\link{pop_init}}
#' @param control.p NULL or named list of buckthorn control treatment parameters 
#'   including:
#'   \describe{
#'     \item{\code{nTrt_grd}}{Proportion of cells with ground treatments in 
#'       each time step}
#'     \item{\code{nTrt_man}}{Proportion of cells with manual treatments in 
#'       each time step}
#'     \item{\code{grd.trt}}{Named vector with ground treatments and associated 
#'       seedling establishment probabilities}
#'     \item{\code{man.trt}}{Named vector with manual treatments and associated
#'       mortality (=success) rates}
#'     \item{\code{add.owners}}{\code{Logical} Do owners treat every year once 
#'       starting a particular treatment?}
#'     \item{\code{t.trt}}{Year to start treatments}
#'   }
#' @return Array of abundances for each cell and age group
#' @keywords run, simulate
#' @export

run_sim <- function(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, control.p=NULL) {
  
  require(tidyverse); require(magrittr)
  
  # Unpack parameters
  tmax <- g.p$tmax
  n.lc <- g.p$n.lc
  dem.st <- g.p$dem.st
  sdd.st <- g.p$sdd.st
  bank <- g.p$bank
  K <- g.p$K  # carrying capacity
  pr.s <- g.p$pr.s  # pre-adult survival
  pr.f <- g.p$pr.f  # pr(fruit)
  fec <- g.p$fec  # mean(fruit per adult)
  age.f <- g.p$age.f  # mean age at first fruiting
  pr.est <- g.p$pr.est  # pr(seedling est)
  pr.sb <- g.p$pr.sb  # pr(ann.surv seed bank)
  lambda <- g.p$lambda  # pop growth rate
  sdd.rate <- g.p$sdd.rate  # 1/mn for dispersal kernel
  pr.eat <- g.p$pr.eat  # pr(birds eat frt)
  pr.s.bird <- g.p$pr.s.bird  # pr(viable | digestion)
  n.ldd <- g.p$n.ldd   # num long distance dispersal events per year
  y.ad <- max(g.p$age.f)
  age.f.d <- length(age.f) > 1
  id_i <- lc.df %>% select(id, id_inbd)
  
  # If buckthorn is being actively managed...
  pr.est.trt <- NULL
  if(!is.null(control.p)) {
    nTrt_grd <- control.p$nTrt_grd * ncell
    nTrt_man <- control.p$nTrt_man * ncell
    grd.trt <- control.p$grd.trt
    man.trt <- control.p$man.trt
    t.trt <- control.p$t.trt
    add.owners <- control.p$add.owners
    est.trt <- tibble(id=numeric(), Trt=character())
    N.trt <- tibble(id=numeric(), Trt=character())
  }
  
  # 1. Initialize populations
  N.sb <- matrix(0, nrow=ngrid, ncol=tmax+1)
  if(age.f.d) {
    N <- array(0, dim=c(ngrid, tmax+1, n.lc, y.ad))  
    N[,1,,] <- N.init
    for(l in 1:n.lc) {
      if(age.f[l] < y.ad) {
        N[,,l,age.f[l]:(y.ad-1)] <- NA
      }
    }
  } else {
    N <- array(0, dim=c(ngrid, tmax+1, y.ad))  
    N[,1,] <- N.init
  }
  
  for(t in 1:tmax) {
    if(age.f.d) {
      N.t <- N[,t,,]
    } else {
      N.t <- N[,t,]
    }
    
    # 2. Implement management
    if(!is.null(control.p) && t >= t.trt) {
      # run functions to implement management controls
      # For future complexity, manual.trt could also push a proportion of
      # the adults back to a previous age so they don't fruit for a number
      # of years after the treatment rather than being killed explicitly
      
      # A. Adjust LC %
      lc.trt <- NULL # to do
      
      # B. Adjust p.est
      if(nTrt_grd > 0) {
        est.trt <- gbPopMod::trt_assign(id_i, ncell, nTrt_grd, grd.trt, 
                                        addOwners=add.owners, trt.m1=est.trt)
        pr.est.trt <- gbPopMod::trt_ground(est.trt, grd.trt)
      }
      
      # C. Adjust N
      if(nTrt_man > 0) {
        N.trt <- gbPopMod::trt_assign(id_i, ncell, nTrt_man, man.trt, 
                                      addOwners=add.owners, trt.m1=N.trt)
        if(age.f.d) {
          N[,t,,] <- gbPopMod::trt_manual(N.t, y.ad, N.trt, man.trt)
        } else {
          N[,t,] <- gbPopMod::trt_manual(N.t, y.ad, N.trt, man.trt)
        }
      }
    }
    
    # 3. Pre-multiply compositional parameters
    pm <- gbPopMod::cell_agg(lc.df, K, pr.s, fec, pr.f, 
                                 pr.eat, pr.est, pr.est.trt)
    
    # 4. Local fruit production
    cat("Year", t, "- Fruiting...")
    N.f <- gbPopMod::make_fruits(N.t, pm$lc.mx, pm$fec.ag, pm$pr.f.ag,
                                  y.ad, age.f.d, dem.st)
    
    # 5. Short distance dispersal
    cat("Dispersing locally...")
    N.seed <- gbPopMod::sdd_disperse(id_i, N.f, pm$pr.eat.ag, pr.s.bird, 
                                sdd.pr, sdd.rate, sdd.st)
    
    # 6. Long distance dispersal
    cat("Dispersing regionally...")
    N.seed <- gbPopMod::ldd_disperse(ncell, id_i, N.seed, n.ldd)
    
    # 7. Seedling establishment
    cat("Establishing...")
    estab.out <- gbPopMod::new_seedlings(ngrid, N.seed, N.sb[,t], pm$pr.est.ag, 
                                          pr.sb, dem.st, bank)
    N.sb[,t+1] <- estab.out$N.sb
    
    # 8. Update abundances
    cat("Updating abundances.\n")
    if(age.f.d) {
      for(l in 1:n.lc) {
        N[,t+1,l,y.ad] <- pmin(round(N[,t,l,y.ad] + 
                                       N[,t,l,age.f[l]-1] * pr.s[l]),
                               pm$K.lc[,l])
        N[,t+1,l,2:(age.f[l]-1)] <- round(N[,t,l,1:(age.f[l]-2)] * pr.s[l])
        N[,t+1,l,1] <- round(estab.out$N.rcrt * pm$rel.dens[,l])
      }
    } else {
      N[,t+1,y.ad] <- pmin(round(N[,t,y.ad] + N[,t,y.ad-1] * pm$pr.s.ag),
                           pm$K.ag)
      N[,t+1,2:(y.ad-1)] <- round(N[,t,1:(y.ad-2)] * pm$pr.s.ag)
      N[,t+1,1] <- estab.out$N.rcrt
    }
  }
  if(age.f.d) {
    N <- apply(N, c(1,2,4), sum, na.rm=TRUE)
  }
  return(list(N=N, N.sb=N.sb))
}
