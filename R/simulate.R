#' Run simulation
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
#' @return Array of abundances for each cell and age group
#' @keywords run, simulate
#' @export

run_sim <- function(ngrid, ncell, g.p, lc.df, sdd.pr, N.init, control.p) {
  
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
  id.i <- lc.df %>% select(id, id.inbd)
  
  # If buckthorn is being actively managed...
  pr.est.trt <- NULL
  if(!is.null(control.p)) {
    nTrt.grd <- control.p$nTrt.grd * ncell
    nTrt.man <- control.p$nTrt.man * ncell
    grd.trt <- control.p$grd.trt
    man.trt <- control.p$man.trt
    lc.chg <- control.p$lc.chg
    n.chg <- control.p$n.chg
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
    cat("Year", t, "")
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
      
      # 2A. Adjust LC %
      if(lc.chg) {
        cat("Changing LC...")
        # i. decide which cells change and how much of each kind of forest
        chg.asn <- cut_assign(n.chg, ncell, lc.df, f.c=6:9)
        
        # ii. cut forest & update SDD neighborhoods
        lc.df[chg.asn$id.chg$id,] <- cut_forest(chg.asn$id.chg, chg.asn$mx, 
                                                f.c=6:9, lc.df)
        sdd.pr[,,,chg.asn$id.chg$id.inbd] <- sdd_set_probs(nrow(chg.asn$id.chg),
                                                           lc.df, g.p, 
                                                           chg.asn$id.chg)
      }
      
      # 2B. Adjust p.est
      if(nTrt.grd > 0) {
        cat("Covering...")
        est.trt <- trt_assign(id.i, ncell, nTrt.grd, grd.trt, 
                                        addOwners=add.owners, trt.m1=est.trt)
        pr.est.trt <- trt_ground(est.trt, grd.trt)
      }
      
      # 2C. Adjust N
      if(nTrt.man > 0) {
        cat("Cutting & spraying...")
        N.trt <- trt_assign(id.i, ncell, nTrt.man, man.trt, 
                                      addOwners=add.owners, trt.m1=N.trt)
        if(age.f.d) {
          N[,t,,] <- trt_manual(N.t, y.ad, N.trt, man.trt)
        } else {
          N[,t,] <- trt_manual(N.t, y.ad, N.trt, man.trt)
        }
      }
    }
    
    # 3. Pre-multiply compositional parameters
    pm <- cell_agg(lc.df, K, pr.s, fec, pr.f, 
                                 pr.eat, pr.est, pr.est.trt)
    
    # 4. Local fruit production
    cat("Fruiting...")
    N.f <- make_fruits(N.t, pm$lc.mx, pm$fec.ag, pm$pr.f.ag,
                                  y.ad, age.f.d, dem.st)
    
    # 5. Short distance dispersal
    cat("Dispersing near...")
    N.seed <- sdd_disperse(id.i, N.f, pm$pr.eat.ag, pr.s.bird, 
                                sdd.pr, sdd.rate, sdd.st)
    
    # 6. Seedling establishment
    cat("Establishing...")
    estab.out <- new_seedlings(ngrid, N.seed, N.sb[,t], pm$pr.est.ag, 
                                          pr.sb, dem.st, bank)
    N.sb[,t+1] <- estab.out$N.sb
    
    # 7. Long distance dispersal
    cat("Dispersing far...")
    estab.out$N.rcrt <- ldd_disperse(ncell, id.i, estab.out$N.rcrt, n.ldd)
    
    # 8. Update abundances
    cat("Updating N...\n")
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
