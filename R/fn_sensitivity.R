#' Run sensitivity analysis
#'
#' This is a wrapper for running simulations over a set of varying parameters.
#' It runs in parallel using the \code{snow} and \code{foreach} packages.
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
#' @param N.init Matrix or array with initial population sizes created by 
#'   \code{\link{pop_init}}
#' @param verbose \code{FALSE} Give updates for each year & process? 
#' @return None
#' @keywords parameters, store, save, output
#' @export

run_sensitivity <- function(p, p.seq, n.sim, ngrid, ncell, g.p, control.p, 
                            lc.df, N.init, verbose=FALSE) {
  require(tidyverse); require(doSNOW); require(foreach)
  
  # set directories
  cat("\n\nSetting directories\n")
  sim.wd <- paste0("out/", ncell, "_t", g.p$tmax, "/")
  par.wd <- paste0(sim.wd, p, "/")
  parSet.wd <- paste0(par.wd, p.seq, "/")
  if(!dir.exists(here("out"))) dir.create(here("out"))
  if(!dir.exists(here(sim.wd))) dir.create(here(sim.wd))
  if(!dir.exists(here(par.wd))) dir.create(here(par.wd))
  
  # dispersal
  cat("\n\nCalculating dispersal neighborhoods\n")
  sdd.pr <- sdd_set_probs(ncell, lc.df, g.p)
  saveRDS(sdd.pr, paste0(par.wd, "sdd_pr.rds"))
  
  cat("\n\nRunning simulations\n")
  sdd.pr <- readRDS(paste0(par.wd, "sdd_pr.rds"))
  for(j in 1:length(p.seq)) {
    # setup for particular parameter value
    set.seed(225)
    g.p[[p]] <- p.seq[[j]]
    
    # run n.sim simulations
    p.c <- makeCluster(g.p$n.cores); registerDoSNOW(p.c)
    out.j <- foreach(i=1:n.sim) %dopar% {
      gbPopMod::run_sim(ngrid, ncell, g.p, lc.df, sdd.pr, 
                        N.init, control.p, verbose=F)
    }
    stopCluster(p.c)
    
    # rearrange output
    ad.j <- map(out.j, ~.$N[,,max(g.p$age.f)]) %>% 
      unlist %>% array(., dim=c(ngrid, g.p$tmax+1, n.sim))
    sb.j <- map(out.j, ~.$N.sb) %>% 
      unlist %>% array(., dim=c(ngrid, g.p$tmax+1, n.sim))
    
    # store output
    if(!dir.exists(here(parSet.wd[j]))) dir.create(here(parSet.wd[j]))
    save_pars(parSet.wd[j], g.p, control.p)
    save_abundances(parSet.wd[j], ad.j, sb.j)
    rm(ad.j); rm(sb.j)
    
    # progress
    cat("Finished parameter set", j, "of", length(p.seq), "\n\n")
  }
  
  
  cat("\nProcessing output\n")
  p.c <- makeCluster(g.p$n.cores); registerDoSNOW(p.c)
  foreach(j=1:length(p.seq), .combine=rbind) %dopar% {
    require(tidyverse); require(stringr); require(gbPopMod)
    options(bitmapType='cairo')
    # setup
    p.j <- paste0(p, ": ", p.seq[j])
    ad.j <- readRDS(paste0(parSet.wd[j], "abund_ad.rds"))
    sb.j <- readRDS(paste0(parSet.wd[j], "abund_sb.rds"))
    g.p <- readRDS(paste0(parSet.wd[j], "pars_glbl.rds"))
    
    # munge data
    ## cell summaries
    ad.mn <- apply(ad.j, 1:2, mean)
    cell.j <- cbind(lc.df, ad.mn) %>% as.tibble %>%
      gather(year, N.adult, (1:ncol(ad.mn)) + ncol(lc.df)) %>%
      rename(ad_Ab=N.adult) %>%
      mutate(ad_sd=c(apply(ad.j, 1:2, sd)),
             ad_pP=c(apply(ad.j > 0, 1:2, mean)),
             sb_Ab=c(apply(log(sb.j+1), 1:2, mean)),
             sb_sd=c(apply(log(sb.j+1), 1:2, sd)),
             sb_pP=c(apply(sb.j > 0, 1:2, mean)),
             ad_L5=c(ad.mn > 0))
    cell.j$ad_L5 <- cell.j$ad_L5 + c(ad.mn > 5)
    cell.j$year <- str_pad(cell.j$year, 3, "left", "0")
    cell.j$p <- p
    if(length(p.seq[[j]])==1) {
      cell.j$p.j <- p.seq[j]
    } else {
      cell.j$p.j <- as.character(p.seq[j])
      cell.j$p.j.OpI <- p.seq[[j]][1]
      cell.j$p.j.Oth <- p.seq[[j]][2]
      cell.j$p.j.Dec <- p.seq[[j]][3]
      cell.j$p.j.WP <- p.seq[[j]][4]
      cell.j$p.j.Evg <- p.seq[[j]][5]
      cell.j$p.j.Mxd <- p.seq[[j]][6]
    }
    ## landscape summaries
    occ.s.ad <- apply(ad.j[lc.df$inbd,,]>0, 2:3, mean)*100
    occ.s.sb <- apply(sb.j[lc.df$inbd,,]>0, 2:3, mean)*100
    less5.s <- apply(ad.j[lc.df$inbd,,]>0 & ad.j[lc.df$inbd,,]<=5, 
                     2:3, mean)*100
    K.ag <- round(as.matrix(lc.df[,4:9]) %*% g.p$K)
    K.s <- apply(ad.j[lc.df$inbd,,]==K.ag[lc.df$inbd,], 
                 2:3, mean)/occ.s.ad*10000
    grid.j <- tibble(year=unique(cell.j$year),
                     pOcc_ad_mn=apply(occ.s.ad, 1, mean),
                     pOcc_ad_sd=apply(occ.s.ad, 1, sd),
                     pOcc_sb_mn=apply(occ.s.sb, 1, mean),
                     pOcc_sb_sd=apply(occ.s.sb, 1, sd),
                     pL5_mn=apply(less5.s, 1, mean),
                     pL5_sd=apply(less5.s, 1, sd),
                     pK_Occ_mn=apply(K.s, 1, mean),
                     pK_Occ_sd=apply(K.s, 1, sd),
                     p=p)
    if(length(p.seq[[j]])==1) {
      grid.j$p.j <- p.seq[j]
    } else {
      grid.j$p.j <- as.character(p.seq[j])
      grid.j$p.j.OpI <- p.seq[[j]][1]
      grid.j$p.j.Oth <- p.seq[[j]][2]
      grid.j$p.j.Dec <- p.seq[[j]][3]
      grid.j$p.j.WP <- p.seq[[j]][4]
      grid.j$p.j.Evg <- p.seq[[j]][5]
      grid.j$p.j.Mxd <- p.seq[[j]][6]
    }
    
    
    # save data summaries
    write_csv(cell.j, paste0(parSet.wd[j], "cell_j.csv"))
    write_csv(grid.j, paste0(parSet.wd[j], "grid_j.csv"))
    rm(ad.j); rm(sb.j)
    
    # save plots
    make_plots_final_t(parSet.wd[j], g.p, filter(cell.j, year==max(cell.j$year)),
                       p.j, 8, 6)
    make_plots_gifs(parSet.wd[j], g.p, cell.j, p.j)
    if(j==1) make_plots_lc(sim.wd, lc.df)
    paste("Finished parameter set", j, "of", length(p.seq))
  }
  stopCluster(p.c)
  
  # grid summary plots
  grid.sum <- map_df(parSet.wd, ~(read_csv(paste0(., "grid_j.csv")))) %>%
    mutate(year=as.numeric(year)) 
  if(length(p.seq[[j]])==1) {
    grid.sum %<>% mutate(p.j=as.factor(p.j)) 
  } else {
    grid.sum %<>% mutate_at(vars(p.j:p.j.Mxd), as.factor)
  }
  make_plots_gridSummary(par.wd, grid.sum, byLC=(length(p.seq[[1]])>1))
}
