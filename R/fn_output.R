#' Save global and control parameters as .rds files
#'
#' This function saves both the global and buckthorn control parameters used in
#' a particular simulation run as a \code{.rds} file in the same directory as
#' the associated plots.
#' @param p.wd Directory to store file in
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @param verbose \code{TRUE} Confirm file creation?
#' @return None
#' @keywords parameters, store, save, output
#' @export

save_pars <- function(p.wd, g.p, control.p, verbose=TRUE) {
  gp_f <- paste0(p.wd, "pars_glbl.rds")
  cp_f <- paste0(p.wd, "pars_ctrl.rds")
  saveRDS(g.p, file=gp_f)
  saveRDS(control.p, file=cp_f)
  if(verbose) {
    if(file.exists(gp_f) & file.exists(cp_f)) { 
      cat("Parameters stored in", p.wd, "\n") 
    } else {
      cat("--- Error: Parameters not stored \n")
    }
  }
}




#' Save simulation output as .rds files
#'
#' This function saves both the global and buckthorn control parameters used in
#' a particular simulation run as a \code{.rds} file in the same directory as
#' the associated plots.
#' @param p.wd Directory to store file in
#' @param ad.N Array of adult abundances with \code{dim=c(ngrid, tmax+1, n.sim)}
#'   generated with a call or calls to \code{\link{run_sim}}
#' @param ad.N Array of seed bank abundances with \code{dim=c(ngrid, tmax+1,
#'   n.sim)} generated with a call or calls to \code{\link{run_sim}}
#' @param verbose \code{TRUE} Confirm file creation?
#' @return None
#' @keywords abundances, store, save, output
#' @export

save_abundances <- function(p.wd, ad.N, sb.N, verbose=TRUE) {
  ad_f <- paste0(p.wd, "abund_ad.rds")
  sb_f <- paste0(p.wd, "abund_sb.rds")
  saveRDS(ad.N, file=ad_f)
  saveRDS(sb.N, file=sb_f)
  if(verbose) {
    if(file.exists(ad_f) & file.exists(sb_f)) { 
      cat("Abundances stored in", p.wd, "\n") 
    } else {
      cat("--- Error: Abundances not stored \n")
    }
  }
}




#' Save plots of final buckthorn distribution
#'
#' This function saves plots of the final distribution of buckthorn. This
#' includes plots for presence/absence and abundance for both seeds (log
#' abundance) and adults.
#' @param p.wd Directory to store file in
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param N.final Dataframe or tibble with processed output from
#'   \code{\link{run_sim}}, including all columns from \code{lc.df} in addition
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.less5}, filtered
#'   to include only the final time step
#' @param txt Text to append to start of plot title
#' @return None
#' @keywords plots, store, save, output
#' @export

make_plots_final_t <- function(p.wd, g.p, N.final, txt=NULL) {
  require(ggplot2); theme_set(theme_bw())
  
  # filenames
  f.ls <- list(paste0(p.wd, "Final_ad_Ab.jpg"), 
               paste0(p.wd, "Final_sb_Ab.jpg"), 
               paste0(p.wd, "Final_ad_sd.jpg"), 
               paste0(p.wd, "Final_sb_sd.jpg"),
               paste0(p.wd, "Final_ad_PA.jpg"),
               paste0(p.wd, "Final_sb_PA.jpg"))
  
  # Adult abundance
  ad.ab.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.adult, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult abundance. Year", g.p$tmax+1))
  ggsave(f.ls[[1]], ad.ab.fin, width=8, height=6, units="in")
  rm(ad.ab.fin)
  
  # Seed bank log abundance
  sb.ab.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.sb, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "log seed abundance. Year", g.p$tmax+1))
  ggsave(f.ls[[2]], sb.ab.fin, width=8, height=6, units="in")
  rm(sb.ab.fin)
  
  # Adult variability
  ad.sd.fin <- ggplot(N.final, aes(x=x, y=-y, fill=sd.ad, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult abundance sd. Year", g.p$tmax+1))
  ggsave(f.ls[[3]], ad.sd.fin, width=8, height=6, units="in")
  rm(ad.sd.fin)
  
  # Seed bank variability
  sb.sd.fin <- ggplot(N.final, aes(x=x, y=-y, fill=sd.sb, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Seed bank abundance sd. Year", g.p$tmax+1))
  ggsave(f.ls[[4]], sb.sd.fin, width=8, height=6, units="in")
  rm(sb.sd.fin)
  
  # Adult presence/absence
  adult.pa.fin <- ggplot(N.final, aes(x=x, y=-y, fill=pP.ad, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult presence. Year", g.p$tmax+1))
  ggsave(f.ls[[5]], adult.pa.fin, width=8, height=6, units="in")
  rm(adult.pa.fin)
  
  # Seed bank presence/absence
  sb.pa.fin <- ggplot(N.final, aes(x=x, y=-y, fill=pP.sb, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Seed presence. Year", g.p$tmax+1))
  ggsave(f.ls[[6]], sb.pa.fin, width=8, height=6, units="in")
  rm(sb.pa.fin)
  
  # Check for success
  for(f in 1:length(f.ls)) {
    if(file.exists(f.ls[[f]])) { 
      cat(f.ls[[f]], "saved\n")
    } else { 
      cat("--- Error:", f.ls[[f]], "not saved! \n") }
  }
}




#' Save gifs of buckthorn distribution through time
#'
#' This function saves gifs of the distribution of buckthorn. This includes
#' plots for presence/absence and abundance for both seeds (log abundance) and
#' adults. There is also a gif distinguishing cells with a low density of
#' adults.
#' @param p.wd Directory to store file in
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param N.out Dataframe or tibble with processed output from
#'   \code{\link{run_sim}}, including all columns from \code{lc.df} in addition
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.less5}
#' @param txt Text to append to start of plot title
#' @return None
#' @keywords plots, gif, store, save, output
#' @export

make_plots_gifs <- function(p.wd, g.p, N.out, txt=NULL) {
  require(gganimate); theme_set(theme_bw())
  
  # filenames
  f.ls <- list(paste0(p.wd, "ad_Ab.gif"),
               paste0(p.wd, "sb_Ab.gif"),
               paste0(p.wd, "ad_sd.gif"),
               paste0(p.wd, "sb_sd.gif"),
               paste0(p.wd, "ad_PA.gif"),
               paste0(p.wd, "sb_PA.gif"),
               paste0(p.wd, "ad_LoDens.gif"))
  
  # Adult abundance
  ad.ab <- ggplot(N.out, aes(x=x, y=-y, fill=N.adult, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult abundance. Year"))
  gganimate(ad.ab, f.ls[[1]], interval=0.2, ani.width=800, ani.height=600)
  rm(ad.ab)
  
  # Seed bank log abundance
  sb.ab <- ggplot(N.out, aes(x=x, y=-y, fill=N.sb, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "log seed abundance. Year"))
  gganimate(sb.ab, f.ls[[2]], interval=0.2, ani.width=800, ani.height=600)
  rm(sb.ab)
  
  # Adult variability
  ad.sd <- ggplot(N.out, aes(x=x, y=-y, fill=sd.ad, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult sd. Year"))
  gganimate(ad.sd, f.ls[[3]], interval=0.2, ani.width=800, ani.height=600)
  rm(ad.sd)
  
  # Seed bank variability
  sb.sd <- ggplot(N.out, aes(x=x, y=-y, fill=sd.sb, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Seed bank sd. Year"))
  gganimate(sb.sd, f.ls[[4]], interval=0.2, ani.width=800, ani.height=600)
  rm(sb.sd)
  
  # Adult presence/absence
  adult.pa <- ggplot(N.out, aes(x=x, y=-y, fill=pP.ad, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Adult presence. Year"))
  gganimate(adult.pa, f.ls[[5]], interval=0.2, ani.width=800, ani.height=600)
  rm(adult.pa)
  
  # Seed bank presence/absence
  sb.pa <- ggplot(N.out, aes(x=x, y=-y, fill=pP.sb, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste(txt, "Seed presence. Year"))
  gganimate(sb.pa, f.ls[[6]], interval=0.2, ani.width=800, ani.height=600)
  rm(sb.pa)
  
  # Low adult density
  adult.lo <- ggplot(N.out, aes(x=x, y=-y, fill=N.less5, frame=year)) +
    geom_tile() + ggtitle(paste(txt, "Blue: Adult density ≤ 5. Year")) + 
    scale_fill_gradient2(low="white", mid="blue", high="pink", midpoint=1)
  gganimate(adult.lo, f.ls[[7]], interval=0.2, ani.width=800, ani.height=600)
  rm(adult.lo)
  
  # Check for success
  for(f in 1:length(f.ls)) {
    if(file.exists(f.ls[[f]])) { 
      cat(f.ls[[f]], "saved\n")
    } else { 
      cat("--- Error:", f.ls[[f]], "not saved! \n") }
  }
}




#' Save plots of initial land cover proportions
#'
#' This function saves plots of the initial coverage of each land cover category
#' @param p.wd Directory to store file in
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @return None
#' @keywords plots, land cover, store, save, output
#' @export

make_plots_lc <- function(p.wd, lc.df) {
  require(ggplot2); theme_set(theme_bw())
  
  # filenames
  f.ls <- list(paste0(p.wd, "LC_OpI.jpg"),
               paste0(p.wd, "LC_Oth.jpg"),
               paste0(p.wd, "LC_Dec.jpg"),
               paste0(p.wd, "LC_WP.jpg"),
               paste0(p.wd, "LC_Evg.jpg"),
               paste0(p.wd, "LC_Mxd.jpg"))
  
  p.opi <- ggplot(lc.df, aes(x=x, y=-y, fill=OpI, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="red", limits=c(0,1))
  ggsave(f.ls[[1]], p.opi, width=10, height=7)
  rm(p.opi)
  
  p.oth <- ggplot(lc.df, aes(x=x, y=-y, fill=Oth, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="gray30", limits=c(0,1))
  ggsave(f.ls[[2]], p.oth, width=10, height=7)
  rm(p.oth)
  
  p.dec <- ggplot(lc.df, aes(x=x, y=-y, fill=Dec, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="green3", limits=c(0,1))
  ggsave(f.ls[[3]], p.dec, width=10, height=7)
  rm(p.dec)
  
  p.wp <- ggplot(lc.df, aes(x=x, y=-y, fill=WP, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="orchid", limits=c(0,1))
  ggsave(f.ls[[4]], p.wp, width=10, height=7)
  rm(p.wp)
  
  p.evg <- ggplot(lc.df, aes(x=x, y=-y, fill=Evg, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="darkgreen", limits=c(0,1))
  ggsave(f.ls[[5]], p.evg, width=10, height=7)
  rm(p.evg)
  
  p.mxd <- ggplot(lc.df, aes(x=x, y=-y, fill=Mxd, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="yellowgreen", limits=c(0,1))
  ggsave(f.ls[[6]], p.mxd, width=10, height=7)
  rm(p.mxd)
  
  # Check for success
  for(f in 1:length(f.ls)) {
    if(file.exists(f.ls[[f]])) { 
      cat(f.ls[[f]], "saved\n")
    } else { 
      cat("--- Error:", f.ls[[f]], "not saved! \n") }
  }
}






