#' Run demographic simulation
#'
#' This function runs the simulation. Currently, it runs all time steps, but for
#' the economic model, the structure will need to be slightly adjusted to run a
#' single time step. The initialization is separated from this function for that
#' reason. 
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd.pr Array with sdd probabilities and neighborhoods created by 
#'   \code{\link{sdd_set_probs}}
#' @param N.init Matrix or array with initial population sizes created by 
#'   \code{\link{pop_init}}
#' @param control.p NULL or named list of buckthorn control treatment parameters 
#'   set with \code{\link{set_control_p}}
#' @param verbose \code{TRUE} Give updates for each year & process? 
#' @return Array of abundances for each cell and age group
#' @keywords run, simulate
#' @export

run_sim <- function(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, 
                    control.p, verbose=TRUE) {
  
  library(tidyverse); library(magrittr)
  
  # Unpack parameters
  list2env(g.p, environment())
  y.ad <- max(age.f)
  id.i <- lc.df %>% select(id, id.in)
  age.f.d <- length(age.f) > 1
  
  # If buckthorn is being actively managed...
  p.est.trt <- NULL
  if(!is.null(control.p)) {
    list2env(control.p, environment())
    nTrt.grd <- ceiling(pTrt.grd * ncell)
    nTrt.man <- ceiling(pTrt.man * ncell)
    nChg <- ceiling(pChg * ncell)
    est.trt <- N.trt <- tibble(id=numeric(), Trt=character())
  }
  
  # 1. Initialize populations
  N.sb <- matrix(0, nrow=ngrid, ncol=tmax+1)
  if(age.f.d) {
    N <- array(0, dim=c(ngrid, tmax+1, n.lc, y.ad))  
    N[,1,,] <- N.init
    for(l in 1:n.lc) {
      if(age.f[l] < y.ad) { N[,,l,age.f[l]:(y.ad-1)] <- NA }
    }
  } else {
    N <- array(0, dim=c(ngrid, tmax+1, y.ad))
    N[,1,] <- N.init
  }
  
  for(t in 1:tmax) {  if(verbose) cat("Year", t, "")
    if(age.f.d) { N.t <- N[,t,,] 
    } else { N.t <- N[,t,] }
    
    # 2. Implement management
    if(!is.null(control.p) && t >= t.trt) {
      # 2A. Adjust LC %
      if(lc.chg && nChg >= 1) {  if(verbose) cat("Change LC...")
        # i. decide which cells change and how much of each kind of forest
        chg.asn <- cut_assign(nChg, ncell, chg.i, lc.df, f.c=6:9)
        # ii. cut forest & update SDD neighborhoods
        lc.df[chg.asn$id.chg$id,] <- cut_forest(chg.asn$id.chg, chg.asn$mx, 
                                                f.c=6:9, lc.df)
        sdd.i <- tibble(id.in=unique(
          arrayInd(which(sdd.pr %in% chg.asn$id.chg$id.in), dim(sdd.pr))[,4]), 
          id=id.i$id[match(id.in, id.i$id.in)])
        sdd.pr[,,,sdd.i$id.in] <- sdd_set_probs(nrow(sdd.i), lc.df, g.p, sdd.i)
      }
      
      # 2B. Adjust p.est
      if(nTrt.grd > 0 || !is.null(grd.i)) {  if(verbose) cat("Cover...")
        est.trt <- trt_assign(id.i=id.i, ncell=ncell, assign_i=grd.i, 
                              nTrt=nTrt.grd, trt.eff=grd.trt, 
                              addOwners=add.owners, trt.m1=est.trt)
        p.est.trt <- trt_ground(est.trt, grd.trt)
      }
      
      # 2C. Adjust N
      if(nTrt.man > 0 || !is.null(man.i)) {  if(verbose) cat("Cut & spray...")
        N.trt <- trt_assign(id.i=id.i, ncell=ncell, assign_i=man.i, 
                            nTrt=nTrt.man, trt.eff=man.trt, 
                            addOwners=add.owners, trt.m1=N.trt)
        if(age.f.d) { N[,t,,] <- trt_manual(N.t, y.ad, N.trt, man.trt)
        } else { N[,t,] <- trt_manual(N.t, y.ad, N.trt, man.trt) }
      }
    }
    
    # 3. Pre-multiply compositional parameters for cell expectations
    pm <- cell_E(lc.df, K, s.jv, fec, p.f, p.eat, p.est, p.est.trt, edges,
                 method)
    
    # 4. Local fruit production
    if(verbose) cat("Fruits...")
    N.f <- make_fruits(N.t, pm$lc.mx, pm$fec.E, pm$p.f.E,
                                  y.ad, age.f.d, dem.st)
    
    # 5. Short distance dispersal
    if(verbose) cat("SDD...")
    N.seed <- sdd_disperse(id.i, N.f, pm$p.eat.E, s.bird, 
                                sdd.pr, sdd.rate, sdd.st, edges)
    
    # 6. Seedling establishment
    if(verbose) cat("Establishment...")
    estab.out <- new_seedlings(ngrid, N.seed, N.sb[,t], pm$p.est.E, 
                               s.sb, dem.st, bank)
    N.sb[,t+1] <- estab.out$N.sb
    
    # 7. Long distance dispersal
    if(verbose) cat("LDD...")
    estab.out$N.rcrt <- ldd_disperse(ncell, id.i, estab.out$N.rcrt, n.ldd)
    
    # 8. Update abundances
    if(verbose) cat("Update N...\n")
    if(age.f.d) {
      for(l in 1:n.lc) {
        N[,t+1,l,y.ad] <- pmin(round(N[,t,l,y.ad] + N[,t,l,age.f[l]-1]*s.jv[l]),
                               pm$K.lc[,l])
        N[,t+1,l,2:(age.f[l]-1)] <- round(N[,t,l,1:(age.f[l]-2)]*s.jv[l])
        N[,t+1,l,1] <- round(estab.out$N.rcrt * pm$rel.dens[,l])
      }
    } else {
      N[,t+1,y.ad] <- pmin(round(N[,t,y.ad] + N[,t,y.ad-1]*pm$s.jv.E), pm$K.E)
      N[,t+1,2:(y.ad-1)] <- round(N[,t,1:(y.ad-2)] * pm$s.jv.E)
      N[,t+1,1] <- estab.out$N.rcrt
    }
  }
  if(age.f.d) N <- apply(N, c(1,2,4), sum, na.rm=TRUE)
  return(list(N=N, N.sb=N.sb))
}






#' Run lambda simulation
#'
#' This function runs the simulation. Currently, it runs all time steps, but for
#' the economic model, the structure will need to be slightly adjusted to run a
#' single time step. The initialization is separated from this function for that
#' reason.
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lambda Vector of length \code{n.lc} with lambdas for each land cover
#'   type
#' @param sdd.pr Array with sdd probabilities and neighborhoods created by
#'   \code{\link{sdd_set_probs}}
#' @param N.init Matrix or array with initial population sizes created by
#'   \code{\link{pop_init}}
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param verbose \code{TRUE} Give updates for each year & process?
#' @return Matrix of abundances for each cell and time step
#' @keywords run, simulate, lambda
#' @export

run_sim_lambda <- function(ngrid, ncell, g.p, lambda, sdd.pr,
                           N.init, verbose=F) {
  library(tidyverse); library(magrittr)
  
  # Unpack parameters
  list2env(g.p, environment())
  id.i <- lc.df %>% select(id, id.in)
  
  # 1. Initialize populations
  N <- matrix(0, ngrid, tmax+1)  
  N[,1] <- apply(N.init, 1, sum)
  
  for(t in 1:tmax){
    # 2. Pre-multiply compositional parameters
    K.E <- as.matrix(lc.df[,4:9]) %*% K
    lambda.E <- as.matrix(lc.df[,4:9]) %*% lambda
    
    # 3. Local growth
    if(verbose) cat("Year", t, "- Grow...")
    N.new <- grow_lambda(N[,t], lambda.E, sdd.rate)
    
    # 4. Short distance dispersal
    if(verbose) cat("SDD...")
    N.emig <- sdd_lambda(N.new, id.i, sdd.pr, sdd.rate, K.E, sdd.st)
    
    # 5. Long distance dispersal
    if(verbose) cat("LDD...")
    N.emig$N <- ldd_disperse(ncell, id.i, N.emig$N, n.ldd)
    
    # 6. Update population sizes
    if(verbose) cat("Update N\n")
    N[,t+1] <- N.emig$N
  }
  return(N)
}
