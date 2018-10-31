#' Set parameters for global sensitivity analysis
#'
#' Generate parameter ranges and details for running a global sensitivity
#' analysis.
#' @param pars Names of parameters to perform sensitivity analysis on
#' @param span \code{"total"} Span for parameter ranges, either \code{"total"}
#'   (i.e., all probabilities [0,1], large ranges for all others, identical
#'   ranges across land cover types) or \code{"gb"} (i.e., plausible ranges for
#'   glossy buckthorn based on experiments, literature, and expert knowledge)
#' @param res \code{"20ac"} Resolution of grid. Either "20ac" or "9km2"
#' @return Named list with a sublist for each parameter, including the parameter
#'   name, the type ('prob', 'cont', 'int'), whether it varies by land cover
#'   category, and the minimum and maximum values
#' @keywords parameters, sensitivity, initialize
#' @export

set_sensitivity_pars <- function(pars, span="total", res="20ac") {
  if(span=="total") {
    par.ls <- list(list(param="p.f", type="prob", LC=1, 
                        min=rep(0,6), max=rep(1,6)),
                   list(param="mu", type="cont", LC=1, 
                        min=rep(0,6), max=rep(500,6)),
                   list(param="gamma", type="cont", LC=0, min=1, max=4),
                   list(param="m", type="int", LC=1, 
                        min=rep(2,6), max=rep(8,6)),
                   list(param="p.c", type="prob", LC=1,
                        min=rep(0,6), max=rep(1,6)),
                   list(param="sdd.rate", type="cont", LC=0, 
                        min=ifelse(res=="20ac", 0.05, 0.5), 
                        max=ifelse(res=="20ac", 0.5, 5)),
                   list(param="sdd.max", type="int", LC=0, 
                        min=ifelse(res=="20ac", 2, 2), 
                        max=ifelse(res=="20ac", 50, 15)),
                   list(param="bird.hab", type="cont", LC=1, 
                        min=rep(0,6), max=rep(1,6)),
                   list(param="n.ldd", type="int", LC=0, min=0, max=10),
                   list(param="s.c", type="prob", LC=0, min=0, max=1),
                   list(param="s.B", type="prob", LC=0, min=0, max=1),
                   list(param="s.M", type="prob", LC=1,
                        min=rep(0,6), max=rep(1,6)),
                   list(param="s.N", type="prob", LC=1,
                        min=rep(0,6), max=rep(1,6)),
                   list(param="K", type="cont", LC=1, 
                        min=rep(0,6), max=rep(40000,6)),
                   list(param="g.D", type="prob", LC=0, min=0, max=1),
                   list(param="g.B", type="prob", LC=0, min=0, max=1),
                   list(param="p", type="prob", LC=1,
                        min=rep(0,6), max=rep(1,6))
    )
  } else if(span=="gb"){
    par.ls <- list(list(param="p.f", type="prob", LC=1, 
                        min=c(0.338, 0, 0.218, 0.232, 0.232, 0.204), 
                        max=c(0.563, 0, 0.364, 0.387, 0.387, 0.340)),
                   list(param="mu", type="cont", LC=1, 
                        min=c(1461, 0, 10, 31, 31, 16), 
                        max=c(2435, 0, 17, 52, 52, 26)),
                   list(param="gamma", type="cont", LC=0, min=2.378, max=2.595),
                   list(param="m", type="int", LC=1, 
                        min=c(2, 2, 4, 4, 4, 4), 
                        max=c(4, 4, 10, 10, 10, 10)),
                   list(param="p.c", type="prob", LC=1,
                        min=c(0.113, 0.113, 0.222, 0.205, 0.205, 0.222), 
                        max=c(0.206, 0.206, 0.314, 0.250, 0.250, 0.314)),
                   list(param="sdd.rate", type="cont", LC=0, 
                        min=ifelse(res=="20ac", 0.107, 1.12), 
                        max=ifelse(res=="20ac", 0.178, 1.87)),
                   list(param="sdd.max", type="int", LC=0, 
                        min=ifelse(res=="20ac", 7, 3), 
                        max=ifelse(res=="20ac", 36, 9)),
                   list(param="bird.hab", type="cont", LC=1, 
                        min=c(0.240, 0.270, 0.037, 0.068, 0.068, 0.068), 
                        max=c(0.400, 0.451, 0.06, 0.113, 0.113, 0.113)),
                   list(param="n.ldd", type="int", LC=0, 
                        min=1, max=ifelse(res=="20ac", 5, 3)),
                   list(param="s.c", type="prob", LC=0, min=0.4875, max=0.605),
                   list(param="s.B", type="prob", LC=0, min=0.641, max=0.766),
                   list(param="s.M", type="prob", LC=1,
                        min=c(0.675, 0, 0.45, 0.45, 0.45, 0.45), 
                        max=c(1, 0.5, 0.75, 0.75, 0.75, 0.75)),
                   list(param="s.N", type="prob", LC=1,
                        min=c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9), 
                        max=c(1, 1, 1, 1, 1, 1)),
                   list(param="K", type="cont", LC=1, 
                        min=c(19533, 0, 1383, 1383, 1383, 1383), 
                        max=c(38713, 100, 5378, 5378, 5378, 5378)),
                   list(param="g.D", type="prob", LC=0, min=0, max=0),
                   list(param="g.B", type="prob", LC=0, min=0.18, max=0.28),
                   list(param="p", type="prob", LC=1,
                        min=c(0.269, 0, 0.316, 0.0594, 0.0594, 0.188), 
                        max=c(0.448, 0, 0.526, 0.099, 0.099, 0.313))
    )
  }
  
  names(par.ls) <- map_chr(par.ls, ~.$param)
  if(res=="9km2") {
    if(span=="total") {
      par.ls$K$min <- rep(0, 6)
      par.ls$K$max <- rep(5000000, 6)
    } else {
      par.ls$K$min <- c(2170287, 0, 153644, 153644, 153644, 153644)
      par.ls$K$max <- c(4301532, 1000, 597578, 597578, 597578, 597578)
    }
  }
  return(par.ls[pars])
}




#' Run global sensitivity analysis
#'
#' Generate and run simulations over a set of varying parameters. It runs in
#' parallel using the \code{snow} and \code{foreach} packages.
#' @param par.ls List output from \code{\link{set_par_ranges}} with parameter to
#'   perform the sensitivity analysis on
#' @param nSamp \code{10} Number of randomly generated parameter sets
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}.
#'   The value of parameter \code{p} will be updated within the loop for each
#'   value in \code{p.seq}.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param sdd \code{NULL} Output with short distance dispersal neighborhoods
#'   created by \code{\link{sdd_set_probs}}. If \code{NULL}, sdd neighborhoods
#'   are calculated
#' @param cell.init Vector with grid id of cells populated at t=0
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param verbose \code{FALSE} Give updates?
#' @param sim.dir \code{"out/sims/"} Directory to store simulation output
#' @return list with elements \code{out} containing the full output from
#'   \link{run_sim}, and \code{results} containing the parameter sets and
#'   simulation summaries
#' @keywords parameters, sensitivity, save, output
#' @export

global_sensitivity <- function(par.ls, nSamp, ngrid, ncell, g.p, lc.df, 
                               sdd=NULL, cell.init, control.p=NULL, 
                               verbose=FALSE, sim.dir="out/sims/") {
  library(tidyverse); library(magrittr); library(foreach); library(doSNOW)
  if(!dir.exists(sim.dir)) dir.create(sim.dir, recursive=TRUE)
  
  # modified from Prowse et al 2016
  if(verbose) cat("Drawing parameters...\n")
  nPar <- length(par.ls)
  samples <- map(g.p[names(par.ls)], 
                 ~matrix(., nrow=nSamp, ncol=length(.), byrow=T))
  raw.samples <- samples
  for(i in 1:nPar) {
    # draw samples from uniform distribution
    raw.samples[[i]][,] <- runif(prod(dim(raw.samples[[i]])), 0, 1)
    # transform to parameter ranges
    samples[[i]][,] <- t(qunif(t(raw.samples[[i]]), 
                               min=par.ls[[i]]$min,
                               max=par.ls[[i]]$max))
    if(par.ls[[i]]$type=="int") {
      samples[[i]][,] <- round(samples[[i]])
    }
  } 
  if(is.null(sdd) && !any(grepl("sdd", names(par.ls)))) {
    sdd <- sdd_set_probs(ncell, lc.df, g.p)
  }
  if(verbose) cat("Running simulations...\n")
  p.c <- makeCluster(g.p$n.cores); registerDoSNOW(p.c)
  out <- foreach(i=1:nSamp, .errorhandling="pass",
                 .packages=c("gbPopMod", "tidyverse", "magrittr")) %dopar% {
    g.p[names(par.ls)] <- map(samples, ~.[i,])
    if(any(grepl("sdd", names(par.ls)))) {
      sdd <- sdd_set_probs(ncell, lc.df, g.p)
    }
    N.init <- pop_init(ngrid, g.p, lc.df, p.0=cell.init, N.0=g.p$N.0)
    sim_i <- run_sim(ngrid, ncell, g.p, lc.df, sdd, N.init, NULL, F, g.p$tmax)
    K.E <- na_if(cell_E(lc.df, g.p$K, g.p$s.M, g.p$s.N, g.p$mu, 
                  g.p$p.f, g.p$p.c, g.p$p)$K.E, 0)
    saveRDS(sim_i$N[,1,dim(sim_i$N)[3]], 
            paste0(sim.dir, "N_", str_pad(i, nchar(nSamp), "left", "0"), ".rds"))
    saveRDS(sim_i$B[,1], 
            paste0(sim.dir, "B_", str_pad(i, nchar(nSamp), "left", "0"), ".rds"))
    saveRDS(sim_i$N[,1,dim(sim_i$N)[3]]/K.E, 
            paste0(sim.dir, "pK_", str_pad(i, nchar(nSamp), "left", "0"), ".rds"))
  }
  stopCluster(p.c)
  
  # calculate grid-wide summaries
  if(verbose) cat("Calculating summaries...\n")
  N <- map(dir(sim.dir, "N_", full.names=T), readRDS)
  B <- map(dir(sim.dir, "B_", full.names=T), readRDS)
  pK <- map(dir(sim.dir, "pK_", full.names=T), readRDS)
  results <- as.data.frame(do.call("cbind", samples))
  par.len <- map_int(samples, ncol)
  par.num <- unlist(list("", paste0("_", 1:6))[(par.len > 1)+1])
  names(results) <- paste0(rep(names(samples), times=par.len), par.num)
  results$pOcc <- map_dbl(N, ~sum(.>0)/ncell)
  results$pSB <- map_dbl(B, ~sum(.>0)/ncell)
  results$pK <- map2_dbl(pK, N, ~mean(.x[.y>0]))
  results$medNg0 <- map_dbl(N, ~median(.[.>0]))
  results$meanNg0 <- map_dbl(N, ~mean(.[.>0]))
  results$sdNg0 <- map_dbl(N, ~sd(.[.>0]))
  return(results)
}




#' Emulate global sensitivity analysis output
#'
#' Emulate the output from a global sensitivity analysis using boosted
#' regression trees with different interaction depths. Based on function
#' described in Prowse et al 2016. Writes files to out/brt, creating /brt if
#' necessary
#' @param sens.out Dataframe of the parameter sets and simulation summaries;
#'   \code{.$results} from \link{global_sensitivity}
#' @param par.ls List output from \code{\link{set_par_ranges}} with parameter to
#'   perform the sensitivity analysis on
#' @param n.cores \code{1} Number of cores for fitting subsample BRTs
#' @param n.sub \code{10} Number of subsamples for each emulation
#' @param td \code{c(1,3,5)} Vector of regression tree interaction depths to
#'   test
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @param brt.dir \code{"out/brt/"} Directory to store boosted regression tree
#'   output
#' @return Success message
#' @keywords parameters, sensitivity, save, output
#' @export

emulate_sensitivity <- function(sens.out, par.ls, n.cores=1, n.sub=10, 
                                td=c(1,3,5), resp, brt.dir="out/brt/") {
  library(tidyverse); library(doSNOW)
  if(!dir.exists(brt.dir)) dir.create(brt.dir)
  x <- which(str_split_fixed(names(sens.out), "_", 2)[,1] %in% names(par.ls))
  y <- which(names(sens.out)==resp)
  sub.prop <- seq(0.5, 1, length.out=n.sub)
  
  p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
  out <- foreach(i=1:n.sub,
                 .packages=c("gbPopMod", "tidyverse")) %dopar% {
    # subset sensitivity results
    sub.samp <- sample_frac(sens.out, sub.prop[i])
    n <- nrow(sub.samp)
    
    # fit BRTs of different tree complexities for given response variable
    for(j in 1:length(td)) {
      td_j <- td[j]
      brt.fit <- dismo::gbm.step(sub.samp, gbm.x=x, gbm.y=y, 
                                 max.trees=200000, n.folds=5, 
                                 family="gaussian", tree.complexity=td_j,
                                 bag.fraction=0.8, silent=T, plot.main=F)
      saveRDS(brt.fit, paste0(brt.dir, resp, "_td-", td_j, "-", n, ".rds"))
    }
  }
  stopCluster(p.c)
  return("Finished fitting BRTs")
}




#' Summarize BRT emulations
#'
#' Summarize the output from global sensitivity analysis boosted regression
#' trees emulations. Based on function described in Prowse et al 2016. Reads BRT
#' output saved via emulate_sensitivity
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @return List of three dataframes: relative influences, cross-validation
#'   deviances, and beta diversity of relative influences (i.e., for stability),
#'   each across different subsample sizes and tree complexities.
#' @param brt.dir \code{"out/brt/"} Directory with stored boosted regression
#'   tree output
#' @keywords parameters, sensitivity, save, output
#' @export

emulation_summary <- function(resp, brt.dir="out/brt/") {
  library(gbm); library(tidyverse)
  f <- dir(brt.dir, paste0(resp, "_"))
  f.i <- str_split_fixed(f, "-", 3)
  cvDev.df <- betaDiv.df <- tibble(response=resp,
                                   td=f.i[,2],
                                   smp=str_remove(f.i[,3], ".rds"))
  cvDev.df$Dev <- NA
  betaDiv.df$beta <- NA
  ri.ls <- vector("list", length(f))
  for(i in seq_along(f)) {
    brt <- readRDS(paste0(brt.dir, f[i]))
    # cross validation deviance
    cvDev.df$Dev[i] <- brt$cv.statistics$deviance.mean
    # relative influence
    ri.ls[[i]] <- as.tibble(brt$contributions) %>%
      mutate(td=cvDev.df$td[i], 
             smp=cvDev.df$smp[i], 
             response=cvDev.df$response[i])
  }
  ri.df <- bind_rows(ri.ls) %>% 
    mutate(param=str_split_fixed(var, "_", 2)[,1]) %>%
    group_by(response, td, smp, param) %>%
    summarise(rel.inf=sum(rel.inf)) %>% 
    ungroup %>% group_by(response, td, smp) %>%
    mutate(rel.inf=rel.inf/sum(rel.inf))
  smp_i <- unique(ri.df$smp)
  for(i in unique(ri.df$td)) {
    for(j in 2:n_distinct(ri.df$smp)) {
      temp <- ri.df %>% ungroup %>% filter(td==i & smp %in% smp_i[c(j-1,j)]) %>%
        spread(smp, rel.inf) %>% dplyr::select(-(1:3)) %>% as.matrix
      beta.div <- as.numeric(MDM::ed(t(temp), q=1, retq=T)['beta'])
      betaDiv.df$beta[betaDiv.df$td==i & betaDiv.df$smp==smp_i[j]] <- beta.div
    }
  }

  return(list(ri.df=ri.df, cvDev.df=cvDev.df, betaDiv.df=betaDiv.df))
}




#' Run univariate sensitivity analyses
#'
#' Running simulations over a set of varying parameters, where each parameter is
#' changed while holding all others constant. It runs in parallel using the
#' \code{doSNOW} and \code{foreach} packages.
#' @param p Parameter to perform the sensitivity analysis on
#' @param p.seq Sequence of values for the parameter
#' @param n.sim \code{2} Number of simulations for each parameter value
#' @param ngrid Number of grid cells in entire map
#' @param ncell Number of inbound grid cells
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}.
#'   The value of parameter \code{p} will be updated within the loop for each
#'   value in \code{p.seq}.
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param verbose \code{FALSE} Give updates for each year & process?
#' @param makeGIFs \code{FALSE} Make a gif for each parameter set?
#' @return Success message
#' @keywords parameters, sensitivity, save, output
#' @export

run_sensitivity <- function(p, p.seq, n.sim, ngrid, ncell, g.p, control.p, 
                            lc.df, verbose=FALSE, makeGIFs=FALSE) {
  library(tidyverse); library(doSNOW); library(foreach)
  
  p.id <- paste(p, collapse="-")
  
  cat("|------------------------------------\n")
  cat("|------ Starting sensitivity analysis for", p.id, "\n")
  cat("|------\n")
  
  # set directories
  byLC <- length(p.seq[[1]])>1
  sim.wd <- paste0("out/", ncell, "_t", g.p$tmax, "/")
  par.wd <- paste0(sim.wd, p.id, "/")
  parSet.wd <- paste0(par.wd, p.seq, "/")
  if(!dir.exists(here("out"))) dir.create(here("out"))
  if(!dir.exists(here(sim.wd))) dir.create(here(sim.wd))
  if(!dir.exists(here(par.wd))) dir.create(here(par.wd))
  
  cat("\nCalculating dispersal neighborhoods\n")
  sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)
  
  cat("\nRunning simulations\n")
  for(j in seq_along(p.seq)) {
    # setup for particular parameter value
    set.seed(225)
    g.p[[p[1]]] <- p.seq[[j]]
    
    # population initialization
    N.init <- pop_init(ngrid, g.p, lc.df)
    
    # dispersal
    if(p[1] %in% c("sdd.max", "sdd.rate", "bird.hab")) {
      sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)
    }
    
    # run n.sim simulations
    p.c <- makeCluster(g.p$n.cores); registerDoSNOW(p.c)
    out.j <- foreach(i=1:n.sim) %dopar% {
      gbPopMod::run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, 
                        N.init, control.p, verbose=F)
    }
    stopCluster(p.c)
    
    # rearrange output
    ad.j <- map(out.j, ~.$N[,,max(g.p$m)]) %>% 
      unlist %>% array(., dim=c(ngrid, g.p$tmax+1, n.sim))
    sb.j <- map(out.j, ~.$B) %>% 
      unlist %>% array(., dim=c(ngrid, g.p$tmax+1, n.sim))
    
    # store output
    if(!dir.exists(here(parSet.wd[j]))) dir.create(here(parSet.wd[j]))
    save_pars(parSet.wd[j], g.p, control.p)
    save_abundances(parSet.wd[j], ad.j, sb.j)
    rm(ad.j); rm(sb.j)
    
    # progress
    cat("  Finished parameter set", j, "of", length(p.seq), "\n\n")
  }
  
  cat("Processing output\n")
  p.c <- makeCluster(g.p$n.cores/3); registerDoSNOW(p.c)
  foreach(j=seq_along(p.seq), .combine=rbind, 
          .packages=c("tidyverse", "stringr", "gbPopMod")) %dopar% {
    options(bitmapType='cairo')
    
    # setup
    p.j <- paste0(p[1], ": ", p.seq[j])
    ad.j <- readRDS(paste0(parSet.wd[j], "abund_ad.rds"))
    sb.j <- readRDS(paste0(parSet.wd[j], "abund_sb.rds"))
    g.p <- readRDS(paste0(parSet.wd[j], "pars_glbl.rds"))
    byLC <- length(p.seq[[j]])>1
    
    # munge data
    ## cell summaries
    ad.mn <- apply(ad.j[lc.df$inbd,,], 1:2, mean)
    K.E <- round(as.matrix(lc.df[,4:9]) %*% g.p$K)
    ### logical: is cell-iter-t at K?
    at.K <- apply(ad.j, c(2,3), function(x) x == K.E & x > 0)
    ### logical: did cell-iter ever reach K?
    reach.K <- at.K
    N.init <- at.K
    for(i in 1:dim(at.K)[1]) {
      for(s in 1:dim(at.K)[3]) {
        reach.K[i,,s] <- sum(reach.K[i,,s]) > 0
      }
      N.init[i,,] <- ad.j[i,1,1] > 0
    }
    ### summarize
    cell_yr.j <- cbind(lc.df[lc.df$inbd,,], ad.mn) %>% as.tibble %>%
      gather(year, N.adult, (1:ncol(ad.mn)) + ncol(lc.df)) %>%
      rename(ad_Ab=N.adult) %>%
      mutate(ad_sd=c(apply(ad.j[lc.df$inbd,,], 1:2, sd)),
             ad_pP=c(apply(ad.j[lc.df$inbd,,] > 0, 1:2, mean)),
             sb_Ab=c(apply(log(sb.j[lc.df$inbd,,]+1), 1:2, mean)),
             sb_sd=c(apply(log(sb.j[lc.df$inbd,,]+1), 1:2, sd)),
             sb_pP=c(apply(sb.j[lc.df$inbd,,] > 0, 1:2, mean)),
             ad_L5=c(ad.mn > 0))
    cell_yr.j$ad_L5 <- cell_yr.j$ad_L5 + c(ad.mn > 5)
    cell_yr.j$year <- str_pad(cell_yr.j$year, 3, "left", "0")
    cell_yr.j$p <- p[1]
    if(byLC) {
      cell_yr.j$p.j <- as.character(p.seq[j])
      cell_yr.j$p.j.Opn <- p.seq[[j]][1]
      cell_yr.j$p.j.Oth <- p.seq[[j]][2]
      cell_yr.j$p.j.Dec <- p.seq[[j]][3]
      cell_yr.j$p.j.WP <- p.seq[[j]][4]
      cell_yr.j$p.j.Evg <- p.seq[[j]][5]
      cell_yr.j$p.j.Mxd <- p.seq[[j]][6]
    } else {
      cell_yr.j$p.j <- p.seq[j]
    }
    ## landscape summaries
    occ.s.ad <- apply(ad.j[lc.df$inbd,,]>0, 2:3, mean)*100
    occ.s.sb <- apply(sb.j[lc.df$inbd,,]>0, 2:3, mean)*100
    less5.s <- apply(ad.j[lc.df$inbd,,]>0 & ad.j[lc.df$inbd,,]<=5, 
                     2:3, mean)*100
    K.s <- apply(at.K[lc.df$inbd,,], 2:3, mean)*100
    t.0K <- apply(ad.j[lc.df$inbd,,]>0 & 
                    !at.K[lc.df$inbd,,] & 
                    reach.K[lc.df$inbd,,] &
                    !N.init[lc.df$inbd,,],
                  c(1,3), sum)
    t.L5 <- apply(ad.j[lc.df$inbd,,]>0 & ad.j[lc.df$inbd,,]<=5,
                  c(1,3), sum)
    grid.j <- tibble(year=unique(cell_yr.j$year),
                     pOcc_ad_mn=apply(occ.s.ad, 1, mean),
                     pOcc_ad_sd=apply(occ.s.ad, 1, sd),
                     pOcc_sb_mn=apply(occ.s.sb, 1, mean),
                     pOcc_sb_sd=apply(occ.s.sb, 1, sd),
                     pL5_mn=apply(less5.s, 1, mean),
                     pL5_sd=apply(less5.s, 1, sd),
                     pK_mn=apply(K.s, 1, mean),
                     pK_sd=apply(K.s, 1, sd),
                     pK_Occ_mn=apply(K.s/occ.s.ad*100, 1, mean),
                     pK_Occ_sd=apply(K.s/occ.s.ad*100, 1, sd),
                     p=p[1])
    cell.j <- cbind(lc.df[lc.df$inbd,], 
                    t0K_mn=apply(t.0K, 1, mean) %>% na_if(0),
                    t0K_sd=apply(t.0K, 1, sd) %>% na_if(0),
                    tL5_mn=apply(t.L5, 1, mean) %>% na_if(0),
                    tL5_sd=apply(t.L5, 1, sd) %>% na_if(0)) %>% as.tibble
    cell.j$p <- p[1]
    if(byLC) {
      grid.j$p.j <- as.character(p.seq[j])
      grid.j$p.j.Opn <- p.seq[[j]][1]
      grid.j$p.j.Oth <- p.seq[[j]][2]
      grid.j$p.j.Dec <- p.seq[[j]][3]
      grid.j$p.j.WP <- p.seq[[j]][4]
      grid.j$p.j.Evg <- p.seq[[j]][5]
      grid.j$p.j.Mxd <- p.seq[[j]][6]
      cell.j$p.j <- as.character(p.seq[j])
      cell.j$p.j.Opn <- p.seq[[j]][1]
      cell.j$p.j.Oth <- p.seq[[j]][2]
      cell.j$p.j.Dec <- p.seq[[j]][3]
      cell.j$p.j.WP <- p.seq[[j]][4]
      cell.j$p.j.Evg <- p.seq[[j]][5]
      cell.j$p.j.Mxd <- p.seq[[j]][6]
    } else {
      grid.j$p.j <- p.seq[j]
      cell.j$p.j <- p.seq[j]
    }
    
    # save data summaries
    write_csv(cell_yr.j, paste0(parSet.wd[j], "cell_yr_j.csv"))
    write_csv(cell.j, paste0(parSet.wd[j], "cell_j.csv"))
    write_csv(grid.j, paste0(parSet.wd[j], "grid_j.csv"))
    rm(ad.j); rm(sb.j)
    
    # save plots
    make_plots_final_t(parSet.wd[j], g.p, 
                       filter(cell_yr.j, year==max(cell_yr.j$year)), p.j, 8, 6)
    if(makeGIFs) make_plots_gifs(parSet.wd[j], g.p, cell_yr.j, p.j)
    if(j==1) make_plots_lc(sim.wd, lc.df)
    paste("  Finished parameter set", j, "of", length(p.seq))
  }
  stopCluster(p.c)
  
  # grid summary plots
  grid.sum <- map_df(parSet.wd, 
                     ~suppressMessages(read_csv(paste0(., "grid_j.csv")))) %>%
    mutate(year=as.numeric(year)) 
  cell.sum <- map_df(parSet.wd, 
                     ~suppressMessages(read_csv(paste0(., "cell_j.csv"))))
  if(byLC) {
    grid.sum %<>% mutate_at(vars(p.j:p.j.Mxd), as.factor)
    cell.sum %<>% mutate_at(vars(p.j:p.j.Mxd), as.factor)
  } else {
    grid.sum %<>% mutate(p.j=as.factor(p.j)) 
    cell.sum %<>% mutate(p.j=as.factor(p.j)) 
  }
  make_plots_gridSummary(p, par.wd, grid.sum, cell.sum, byLC)
  
  # progress
  cat("\n")
  cat("|------\n")
  cat("|------ Finished sensitivity analysis for", p.id, "\n")
  cat("|------------------------------------\n\n\n")
  return(paste("Completed", p.id))
}
