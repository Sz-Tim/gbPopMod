#' Run demographic simulation
#'
#' Run the simulation. Currently, it runs all time steps, but for the economic
#' model, the structure will need to be slightly adjusted to run a single time
#' step. The initialization is separated from this function for that reason.
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd Output with short distance dispersal neighborhoods created by
#'   \code{\link{sdd_set_probs}}
#' @param N.init Matrix or array with initial population sizes created by
#'   \code{\link{pop_init}}
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param verbose \code{TRUE} Show progress bar?
#' @param save_yrs \code{NULL} Vector of years to save; if \code{NULL}, all
#'   years are returned
#' @param K_max \code{NULL} Maximum carrying capacity to avoid memory overloads
#' @param dem_out \code{TRUE} Store nFl, nSd, nSdStay, D?
#' @return Array N of abundances for each cell and age group, matrix B of seed
#'   bank abundances, matrix nFl with number of flowering individuals, matrix
#'   nSd of total seeds produced in each cell, matrix nSdStay of total seeds
#'   remaining in each cell, and matrix D of immigrants into each cell.
#' @keywords run, simulate
#' @export

run_sim <- function(ngrid, ncell, g.p, lc.df, sdd, N.init, control.p, 
                    verbose=TRUE, save_yrs=NULL, K_max=NULL, dem_out=TRUE) {
  
  library(tidyverse); library(magrittr)
  
  # Unpack parameters
  list2env(g.p, environment())
  m.max <- max(m)
  id.i <- lc.df %>% select(id, id.in)
  m.d <- n_distinct(m) > 1  # does m vary across cells?
  if(method=="lm" && m.d) stop("m must be uniform if method==\"lm\"")
  
  # identify which years to store data
  if(is.null(save_yrs)) {
    n_yrs <- tmax
    yrs <- setNames(1:tmax, 1:tmax)
  } else {
    n_yrs <- length(save_yrs)
    yrs <- setNames(1:n_yrs, save_yrs)
  }
  
  # If buckthorn is being actively managed...
  p.trt <- NULL
  if(!is.null(control.p)) {
    list2env(control.p, environment())
    est.trt <- N.trt <- tibble(id=numeric(), Trt=character())
  }
  
  # 1. Initialize populations
  ## storage objects to return: x
  ## annual temp objects: x.0 = start of year, x.1 = end of year
  N.0 <- N.1 <- N.init
  B.0 <- rep(0, ngrid)
  B <- matrix(0, nrow=ngrid, ncol=n_yrs)
  if(dem_out) nFl <- nSd <- nSdStay <- D <- B
  if(m.d) {
    N <- array(0, dim=c(ngrid, n_yrs, n.lc, m.max))  
  } else {
    N <- array(0, dim=c(ngrid, n_yrs, m.max))
  }
  # 2. Pre-multiply compositional parameters for cell expectations
  pm <- cell_E(lc.df, K, s.M, s.N, mu, p.f, p.c, p, p.trt, edges, method)
  if(!is.null(K_max)) pm$K.E <- pmin(K_max, pm$K.E)
  
  if(verbose) pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
  for(k in 1:tmax) {  
    
    if(dem_out) nFl.1 <- nSd.1 <- nSdStay.1 <- D.1 <- rep(0, ngrid)
    
    # 3. Implement management
    if(!is.null(control.p) && k >= t.trt) {
      # 3A. Adjust LC %
      if(lc.chg && nChg >= 1) { 
        # i. decide which cells change and how much of each kind of forest
        chg.asn <- cut_assign(pChg, ncell, chg.i, lc.df, forest.col=6:9)
        # ii. cut forest & update SDD neighborhoods
        lc.df[chg.asn$id.chg$id,] <- cut_forest(chg.asn$id.chg, chg.asn$mx, 
                                                forest.col=6:9, lc.df)
        sdd.i <- tibble(id.in=unique(
          arrayInd(which(sdd$i %in% chg.asn$id.chg$id.in), dim(sdd$i))[,4]), 
          id=id.i$id[match(id.in, id.i$id.in)])
        sdd_new <- sdd_update_probs(lc.df, g.p, sdd.i, sdd$i)
        sdd$i[,,,sdd.i$id.in] <- sdd_new$i
        sdd$sp[sdd.i$id.in] <- sdd_new$sp
      }
      
      # 3B. Adjust p
      if(pTrt.grd*ncell > 0.5 || !is.null(grd.i)) { 
        est.trt <- trt_assign(id.i=id.i, ncell=ncell, assign_i=grd.i, 
                              pTrt=pTrt.grd, trt.eff=grd.trt, 
                              addOwners=add.owners, trt.m1=est.trt)
        p.trt <- trt_ground(est.trt, grd.trt)
        pm <- cell_E(lc.df, K, s.M, s.N, mu, p.f, p.c, p, p.trt, edges, method)
        if(!is.null(K_max)) pm$K.E <- pmin(K_max, pm$K.E)
      }
      
      # 3C. Adjust N
      if(pTrt.man*ncell > 0.5 || !is.null(man.i)) { 
        N.trt <- trt_assign(id.i=id.i, ncell=ncell, assign_i=man.i, 
                            pTrt=pTrt.man, trt.eff=man.trt, 
                            addOwners=add.owners, trt.m1=N.trt)
        N.0 <- trt_manual(N.0, m.max, N.trt, man.trt)
      }
    }
    
    # 4. Local fruit production
    N.f <- make_fruits(N.0, pm$lc.mx, pm$mu.E, pm$p.f.E, m.max, m.d, dem.st)
    if(dem_out) nFl.1[N.f$id] <- N.f$N.rpr
    
    # 5. Short distance dispersal
    N.Sd <- sdd_disperse(id.i, N.f, gamma, pm$p.c.E, s.c, 
                         sdd$sp, sdd.rate, sdd.st, edges)
    rm(N.f)
    if(dem_out) {
      nSd.1[N.Sd$N.source$id] <- N.Sd$N.source$N.produced
      nSdStay.1[N.Sd$N.source$id] <- N.Sd$N.source$N.dep
      D.1[N.Sd$N.seed$id] <- N.Sd$N.seed$N
      D.1 <- D.1 - nSdStay.1
    }
    
    # 6. Seedling establishment
    estab.out <- new_seedlings(ngrid, N.Sd$N.seed, B.0, pm$p.E, g.D, g.B,
                               s.B, dem.st, bank)
    rm(N.Sd)
    B.1 <- estab.out$B
    
    # 7. Long distance dispersal
    estab.out$M.0 <- ldd_disperse(ncell, id.i, estab.out$M.0, n.ldd)
    
    # 8. Calculate abundances
    N.1[] <- 0
    if(m.d) {
        for(l in 1:n.lc) {
          N.1[,l,1] <- round(estab.out$M.0 * pm$s.M.E)
          N.1[,l,2:(m[l]-1)] <- round(N.0[,l,1:(m[l]-2)]*s.M[l])
          N.1[,l,m.max] <- pmin(round(N.0[,l,m.max]*s.N[l] + 
                                        N.0[,l,m[l]-1]*s.M[l]),
                                pm$K.lc[,l])
        }
    } else if(m.max > 2) {
      N.1[,1] <- round(estab.out$M.0 * pm$s.M.E)
      N.1[,2:(m.max-1)] <- round(N.0[,1:(m.max-2)]*pm$s.M.E)
      N.1[,m.max] <- pmin(round(N.0[,m.max]*pm$s.N.E + N.0[,m.max-1]*pm$s.M.E), 
                          pm$K.E)
    } else if(m.max == 2 ) {
      N.1[,1] <- round(estab.out$M.0 * pm$s.M.E)
      N.1[,m.max] <- pmin(round(N.0[,m.max]*pm$s.N.E + N.0[,1]*pm$s.M.E), 
                          pm$K.E)
    } else {
      N.1[,1] <- pmin(round(N.0[,1]*pm$s.N.E + estab.out$N.rcrt*pm$s.M.E), 
                      pm$K.E)
    }
    rm(estab.out)
    
    # 9. Store years
    if(k %in% names(yrs)) {
      yr_k <- which(names(yrs)==k)
      if(m.d) { N[,yr_k,,] <- N.1
      } else { N[,yr_k,] <- N.1 }
      B[,yr_k] <- B.1
      if(dem_out) {
        nFl[,yr_k] <- nFl.1
        nSd[,yr_k] <- nSd.1
        nSdStay[,yr_k] <- nSdStay.1
        D[,yr_k] <- D.1
      }
    }
    
    # 10. Update objects
    N.0 <- N.1
    B.0 <- B.1
    
    if(verbose) setTxtProgressBar(pb, k)
  }
  if(verbose) close(pb)
  if(m.d) N <- apply(N, c(1,2,4), sum, na.rm=TRUE)
  if(!dem_out) nFl <- nSd <- nSdStay <- D <- NULL
  return(list(N=N, B=B, nFl=nFl, nSd=nSd, nSdStay=nSdStay, D=D))
}






#' Run lambda simulation
#'
#' Run the simulation. Currently, it runs all time steps, but for the economic
#' model, the structure will need to be slightly adjusted to run a single time
#' step. The initialization is separated from this function for that reason.
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lambda Vector of length \code{n.lc} with lambdas for each land cover
#'   type
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd.pr Array with sdd probabilities and neighborhoods created by
#'   \code{\link{sdd_set_probs}}
#' @param N.init Matrix or array with initial population sizes created by
#'   \code{\link{pop_init}}
#' @param method \code{"wt.mn"} Method for calculating cell expectations, taking
#'   values of \code{"wt.mn"} or \code{"lm"}. If \code{"wt.mn"}, the expectation
#'   for each parameter is the weighted mean across land cover types
#'   proportional to their coverage, with the land cover specific values stored
#'   in the parameter vectors. If \code{"lm"}, the expectation is calculated in
#'   a regression with the slopes contained in each parameter vector.
#'   Individuals cannot be assigned to specific land cover categories with
#'   \code{"lm"}, so \code{"m"} must be scalar.
#' @param verbose \code{TRUE} Show progress bar?
#' @return Matrix N of abundances for each cell and time step and vector
#'   lambda.E of predicted lambda value for each cell
#' @keywords run, simulate, lambda
#' @export

run_sim_lambda <- function(ngrid, ncell, g.p, lambda, lc.df, sdd.pr,
                           N.init, method="wt.mn", verbose=F) {
  library(tidyverse); library(magrittr)
  
  # Unpack parameters
  list2env(g.p, environment())
  id.i <- lc.df %>% select(id, id.in)
  
  # 1. Initialize populations
  N <- matrix(0, ngrid, tmax+1)  
  N[,1] <- apply(N.init, 1, sum)
  
  if(verbose) pb <- txtProgressBar(min=0, max=tmax, width=80, style=3)
  for(t in 1:tmax){
    # 2. Pre-multiply compositional parameters
    if(method=="wt.mn") {
      K.E <- as.matrix(lc.df[,4:9]) %*% K
      lambda.E <- as.matrix(lc.df[,4:9]) %*% lambda
    } else if(method=="lm") {
      lc.mx <- cbind(1, as.matrix(select(lc.df, 
                              -one_of("x", "y", "x_y", "inbd", "id", "id.in"))))
      K.E <- exp(lc.mx[,1:length(K)] %*% K)
      lambda.E <- exp(lc.mx[,1:length(lambda)] %*% lambda)
    } 
    
    # 3. Local growth
    N.new <- grow_lambda(N[,t], lambda.E, sdd.rate)
    
    # 4. Short distance dispersal
    N.emig <- sdd_lambda(N.new, id.i, sdd.pr, sdd.rate, K.E, sdd.st)
    
    # 5. Long distance dispersal
    N.emig$N <- ldd_disperse(ncell, id.i, N.emig$N, n.ldd)
    
    # 6. Update population sizes
    N[,t+1] <- N.emig$N
    if(verbose) setTxtProgressBar(pb, t)
  }
  if(verbose) close(pb)
  return(list(N=N, lam.E=lambda.E))
}




#' Implement management and iterate buckthorn populations for one time step
#'
#' Run one time step of the simulation, possibly implementing management actions
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param N.0 Array with initial population sizes, either returned from a
#'   previous iteration or read from a stored .rds file in \code{path}
#' @param B.0 Vector of initial seed bank abundances, either returned from a
#'   previous iteration or read from a stored .rds file in \code{path}
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd Output with short distance dispersal neighborhoods created by
#'   \code{\link{sdd_set_probs}}
#' @param control.p \code{NULL} NULL or named list of buckthorn control
#'   treatment parameters set with \code{\link{set_control_p}}
#' @param grd_cover.i \code{NULL} Dataframe with a row for each cell
#'   implementing ground cover management, and columns \code{id} and \code{Trt}
#'   detailing the cell ID and ground cover treatment (\code{"Cov", "Com",
#'   "Lit"} for ground cover crop, compaction, or litter, respectively).
#' @param mech_chem.i \code{NULL} Dataframe with a row for each cell
#'   implementing manual management of adults, and columns \code{id} and
#'   \code{Trt} detailing the cell ID and manual treatment (\code{"M", "C",
#'   "MC"} for mechanical, chemical, or mechanical and chemical, respectively).
#' @param read_write \code{FALSE} Read and write \code{N} and \code{B}
#' @param path \code{NULL} Directory for stored output. Overwrites files
#'   (path/N.rds, path/B.rds) each iteration.
#' @return Array N of abundances for each cell and age group, and vector B of
#'   seed bank abundances.
#' @keywords run, simulate
#' @export

iterate_pop <- function(ngrid, ncell, N.0=NULL, B.0=NULL, g.p, lc.df, sdd, 
                        control.p=NULL, grd_cover.i=NULL, mech_chem.i=NULL, 
                        read_write=FALSE, path=NULL) {
  library(gbPopMod); library(tidyverse); library(magrittr)
  if(read_write) {
    N.0 <- readRDS(paste0(path, "/N.rds"))
    B.0 <- readRDS(paste0(path, "/B.rds"))
  }
  list2env(g.p, environment())
  id.i <- lc.df %>% dplyr::select(id, id.in)
  m.max <- max(m)
  N.1 <- array(0, dim=dim(N.0))
  
  #--- update parameters
  if(!is.null(grd_cover.i)) {
    p.trt <- trt_ground(grd_cover.i, control.p$grd.trt)
  } else {
    p.trt <- NULL
  }
  if(!is.null(mech_chem.i)) {
    N.0 <- trt_manual(N.0, m.max, mech_chem.i, control.p$man.trt, 
                      ncol(mech_chem.i)>2)
  }
  # pre-multiply compositional parameters for cell expectations
  pm <- cell_E(lc.df, K, s.M, s.N, mu, p.f, p.c, p, p.trt)
  
  #--- fruit production
  N.f <- make_fruits(N.0, pm$lc.mx, pm$mu.E, pm$p.f.E, m.max, T)
  
  #--- short distance dispersal
  N.Sd <- sdd_disperse(id.i, N.f, gamma, pm$p.c.E, s.c, sdd$sp, sdd.rate)
  
  #--- seedling establishment
  estab.out <- new_seedlings(ngrid, N.Sd$N.seed, B.0, pm$p.E, g.D, g.B, s.B)
  estab.out$M.0 <- ldd_disperse(ncell, id.i, estab.out$M.0, n.ldd)
  
  #--- update abundances
  B.1 <- estab.out$B
  for(l in 1:6) {
    N.1[,l,1] <- ceiling(estab.out$M.0 * pm$lc.mx[,l])
    N.1[,l,2:(m[l]-1)] <- round(N.0[,l,1:(m[l]-2)]*s.M[l])
    N.1[,l,m.max] <- pmin(round(N.0[,l,m.max]*s.N[l] + N.1[,l,m[l]-1]*s.M[l]),
                          pm$K.lc[,l])
  }
  if(read_write) {
    saveRDS(N.1, paste0(path, "/N.rds"))
    saveRDS(B.1, paste0(path, "/B.rds"))
  }
  return(list(N=N.1, B=B.1))
}





#' Implement management and iterate buckthorn for one time step for parcels
#'
#' Run one time step of the simulation, implementing management actions at the
#' sub-pixel parcel level. This is specifically designed for integration with
#' the USDA-NIFA economic decision model, where management actions are taken by
#' individual parcels which may be smaller than the land cover map pixels. The
#' indexing is consequently different than for the other iteration functions.
#' @param parcel.df Dataframe with index columns ("id", "id.in", "id.pp"), and
#'   columns detailing land cover proportions and grid proportions.
#' @param pp.ls List of length ncell where each element i identifies which
#'   id.i$id.pp are within pixel i
#' @param N.0 Array with initial population sizes, either returned from a
#'   previous iteration or read from a stored .rds file in \code{path}
#' @param B.0 Vector of initial seed bank abundances, either returned from a
#'   previous iteration or read from a stored .rds file in \code{path}
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd Output with short distance dispersal neighborhoods created by
#'   \code{\link{sdd_set_probs}}
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param grd_cover.i \code{NULL} Dataframe with a row for each cell
#'   implementing ground cover management, and columns \code{id} and \code{Trt}
#'   detailing the cell ID and ground cover treatment (\code{"Cov", "Com",
#'   "Lit"} for ground cover crop, compaction, or litter, respectively).
#' @param mech_chem.i \code{NULL} Dataframe with a row for each cell
#'   implementing manual management of adults, and columns \code{id} and
#'   \code{Trt} detailing the cell ID and manual treatment (\code{"M", "C",
#'   "MC"} for mechanical, chemical, or mechanical and chemical, respectively).
#' @param read_write \code{FALSE} Read and write \code{N} and \code{B}
#' @param path \code{NULL} Directory for stored output. Overwrites files
#'   (path/N.rds, path/B.rds) each iteration.
#' @return Array N of abundances for each cell and age group, and vector B of
#'   seed bank abundances.
#' @keywords run, simulate
#' @export

iterate_pop_econ <- function(parcel.df, pp.ls, N.0, B.0, g.p, lc.df, sdd, 
                             control.p, grd_cover.i=NULL, mech_chem.i=NULL, 
                             read_write=FALSE, path=NULL) {
  library(gbPopMod); library(tidyverse); library(magrittr)
  if(read_write) {
    N.0 <- readRDS(paste0(path, "/N.rds"))
    B.0 <- readRDS(paste0(path, "/B.rds"))
  }
  LCs <- names(lc.df)[(1:dim(N.0)[2])+3]
  ncell <- n_distinct(parcel.df$id.in)
  npp <- nrow(parcel.df)
  list2env(g.p, environment())
  m.max <- max(m)
  N.1 <- array(0, dim=dim(N.0))
  
  #--- pre-multiply compositional parameters for cell expectations
  p.trt <- NULL
  pm <- cell_E(lc.df, K, s.M, s.N, mu, p.f, p.c, p, p.trt)
  
  #--- implement management
  if(!is.null(grd_cover.i)) {
    p.trt <- trt_ground(grd_cover.i, control.p$grd.trt)
  }
  if(!is.null(mech_chem.i)) {
    N.0 <- trt_manual(N.0, m.max, mech_chem.i, control.p$man.trt, 
                      ncol(mech_chem.i)>2)
  }
  
  #--- sum abundance within pixels
  N.0.px <- unique(parcel.df$id.in[rowSums(N.0[,,7])>0])
  N.px <- array(0, dim=c(ncell, length(LCs), m.max))
  N.px[N.0.px,,] <- aperm(vapply(N.0.px, 
                                 function(x) apply(N.0[pp.ls[[x]],,], 2:3, sum),
                                 N.0[1,,]), c(3,1,2))
  
  #--- fruit production
  N.f <- make_fruits(N.px, pm$lc.mx[lc.df$inbd,], pm$mu.E[lc.df$inbd,,drop=F], 
                     pm$p.f.E[lc.df$inbd,,drop=F], m.max, T)
  N.f$id <- lc.df$id[match(N.f$id, lc.df$id.in)]  # convert from id.in to id
  
  #--- short distance dispersal
  N.Sd <- sdd_disperse(lc.df[,c("id", "id.in")], N.f, gamma, pm$p.c.E, 
                       s.c, sdd$sp, sdd.rate)
  rm(N.f)
  
  #--- seedling establishment
  estab.out <- new_seedlings(ngrid, N.Sd$N.seed, B.0, pm$p.E, g.D, g.B, s.B)
  rm(N.Sd)
  estab.out$M.0 <- estab.out$M.0[lc.df$inbd]
  
  #--- allocate seedlings among parcels
  B.1 <- estab.out$B
  for(l in 1:n.lc) {
    N.1[,l,1] <- round(estab.out$M.0[parcel.df$id.in] * 
                         parcel.df[,LCs[l]] * parcel.df$Grid_Proportion)
    N.1[,l,2:(m[l]-1)] <- round(N.0[,l,1:(m[l]-2)]*s.M[l])
    N.1[,l,m.max] <- pmin(round(N.0[,l,m.max]*s.N[l] + N.1[,l,m[l]-1]*s.M[l]),
                          parcel.df[,paste0("K_", LCs[l])])
  }
  rm(estab.out)
  
  #--- retroactively apply ground cover treatment by recalculating establishment
  if(!is.null(p.trt)) {
    N.1[p.trt$id,,1] <- round(N.1[p.trt$id,,1] / 
                                pm$p.E[parcel.df$id[p.trt$id]] * p.trt$p)
  }
  
  #--- long distance dispersal
  if(n.ldd > 0) {
    ldd_id.pp <- sample.int(npp, n.ldd)
    N.1[ldd_id.pp,, 1] <- N.1[ldd_id.pp,, 1] + 1
  }
  
  if(read_write) {
    saveRDS(N.1, paste0(path, "/N.rds"))
    saveRDS(B.1, paste0(path, "/B.rds"))
  }
  return(list(N=N.1, B=B.1))
}

