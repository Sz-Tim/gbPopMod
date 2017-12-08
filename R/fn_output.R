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
      cat("  Parameters stored in", p.wd, "\n") 
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
      cat("  Abundances stored in", p.wd, "\n") 
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
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.L5}, filtered
#'   to include only the final time step
#' @param txt Text to append to start of plot title
#' @param w Figure output width (inches)
#' @param h Figure output height (inches)
#' @return None
#' @keywords plots, store, save, output
#' @export

make_plots_final_t <- function(p.wd, g.p, N.final, txt=NULL, w=8, h=6) {
  library(ggplot2); theme_set(theme_bw())
  
  # filenames
  f.nm <- c("ad_Ab", "sb_Ab", "ad_sd", "sb_sd", "ad_pP", "sb_pP")
  f.full <- paste0(p.wd, "Final_", f.nm, ".jpg") 
  
  p.fin <- ggplot(N.final, aes(x=x, y=-y, colour=inbd)) +
    scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA))
  
  # Adult abundance
  ggsave(f.full[1], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[1])) +
           ggtitle(paste(txt, "Adult mean(abundance). Year", g.p$tmax+1)))
  ggsave(f.full[2], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[2])) +
           ggtitle(paste(txt, "Seed mean(log abundance). Year", g.p$tmax+1)))
  ggsave(f.full[3], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[3])) +
           ggtitle(paste(txt, "Adult sd(abundance). Year", g.p$tmax+1)))
  ggsave(f.full[4], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[4])) +
           ggtitle(paste(txt, "Seed sd(log abundance). Year", g.p$tmax+1)))
  ggsave(f.full[5], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[5])) +
           ggtitle(paste(txt, "Adult mean(presence). Year", g.p$tmax+1)))
  ggsave(f.full[6], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[6])) +
           ggtitle(paste(txt, "Seed mean(presence). Year", g.p$tmax+1)))
  
  # Check for success
  for(f in 1:length(f.full)) {
    if(file.exists(f.full[f])) { 
      cat(" ", f.full[f], "saved\n")
    } else { 
      cat("--- Error:", f.full[f], "not found! \n") }
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
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.L5}
#' @param txt Text to append to start of plot title
#' @param w Figure output width (pixels)
#' @param h Figure output height (pixels)
#' @param i Gif framerate
#' @return None
#' @keywords plots, gif, store, save, output
#' @export

make_plots_gifs <- function(p.wd, g.p, N.out, txt=NULL, w=800, h=600, i=0.2) {
  library(gganimate); theme_set(theme_bw())
  
  # filenames
  f.nm <- c("ad_Ab", "sb_Ab", "ad_sd", "sb_sd", "ad_pP", "sb_pP", "ad_L5")
  f.full <- paste0(p.wd, f.nm, ".gif")
  
  p.gif <- ggplot(N.out, aes(x=x, y=-y, frame=year, colour=inbd)) + 
    scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA))
  
  # Adult abundance
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[1])) + 
              ggtitle(paste(txt, "Adult mean(abundance). Year")), 
            filename=f.full[1], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[2])) + 
              ggtitle(paste(txt, "Seed mean(log abundance). Year")), 
            filename=f.full[2], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[3])) + 
              ggtitle(paste(txt, "Adult sd(abundance). Year")), 
            filename=f.full[3], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[4])) + 
              ggtitle(paste(txt, "Seed sd(log abundance). Year")), 
            filename=f.full[4], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[5])) + 
              ggtitle(paste(txt, "Adult mean(presence). Year")), 
            filename=f.full[5], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[6])) + 
              ggtitle(paste(txt, "Seed mean(presence). Year")), 
            filename=f.full[6], interval=i, ani.width=w, ani.height=h)
  gganimate(p=p.gif + geom_tile(aes_string(fill=f.nm[7])) %+% 
              scale_fill_gradient2(low="white", mid="blue", 
                                   high="pink", midpoint=1) +
              ggtitle(paste(txt, "Adult mean(0 < abundance ≤ 5). Year")), 
            filename=f.full[7], interval=i, ani.width=w, ani.height=h)
  
  # Check for success
  for(f in 1:length(f.full)) {
    if(file.exists(f.full[f])) { 
      cat(" ", f.full[f], "saved\n")
    } else { 
      cat("--- Error:", f.full[f], "not found! \n") }
  }
}




#' Save plots of initial land cover proportions
#'
#' This function saves plots of the initial coverage of each land cover category
#' @param p.wd Directory to store file in
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @param w Figure output width (inches)
#' @param h Figure output height (inches)
#' @return None
#' @keywords plots, land cover, store, save, output
#' @export

make_plots_lc <- function(p.wd, lc.df, w=10, h=7) {
  library(ggplot2); theme_set(theme_bw())
  
  # filenames
  f.full <- paste0(p.wd, c("LC_OpI.jpg", "LC_Oth.jpg", "LC_Dec.jpg", 
                           "LC_WP.jpg", "LC_Evg.jpg", "LC_Mxd.jpg"))
  
  p.lc <- ggplot(lc.df, aes(x=x, y=-y, colour=inbd)) + 
    scale_colour_manual(values=c("gray", NA))
  
  ggsave(f.full[1], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=OpI)) + ggtitle("Open Invasible") +
           scale_fill_gradient(low="white", high="red", limits=c(0,1)))
  ggsave(f.full[2], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Oth)) + ggtitle("Open Invasible")+
           scale_fill_gradient(low="white", high="gray30", limits=c(0,1)))
  ggsave(f.full[3], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Dec)) + ggtitle("Open Invasible")+
           scale_fill_gradient(low="white", high="green3", limits=c(0,1)))
  ggsave(f.full[4], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=WP)) + ggtitle("Open Invasible")+
           scale_fill_gradient(low="white", high="orchid", limits=c(0,1)))
  ggsave(f.full[5], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Evg)) + ggtitle("Open Invasible")+
           scale_fill_gradient(low="white", high="darkgreen", limits=c(0,1)))
  ggsave(f.full[6], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Mxd)) + ggtitle("Open Invasible")+
           scale_fill_gradient(low="white", high="yellowgreen", limits=c(0,1)))
  
  # Check for success
  for(f in 1:length(f.full)) {
    if(file.exists(f.full[f])) { 
      cat(" ", f.full[f], "saved\n")
    } else { 
      cat("--- Error:", f.full[f], "not found! \n") }
  }
}




#' Save plots of buckthorn distribution summaries
#'
#' This function saves plots of buckthorn distribution summary statistics for
#' the whole grid across simulations. This includes plots for the mean and sd of
#' the proportion of cells occupied, of the proportion of low density cells, and
#' of proportion of occupied cells at carrying capacity.
#' @param par.wd Directory to store file in
#' @param grid.sum Dataframe or tibble with processed output for whole-grid
#'   summaries
#' @param byLC \code{FALSE} For LC-specific parameters, should separate plots be
#'   created for each LC?
#' @param txt Text to append to start of plot title
#' @param w Figure output width (inches)
#' @param h Figure output height (inches)
#' @return None
#' @keywords plots, store, save, output
#' @export

make_plots_gridSummary <- function(p.wd, grid.sum, byLC=FALSE, txt=NULL,
                                   w=8, h=6) {
  library(scales); library(ggplot2); theme_set(theme_bw())
  p <- grid.sum$p[1]
  n.set <- length(unique(grid.sum$p.j))
  p.col <- seq_gradient_pal(low="#e5f5e0", high="#00441b")((1:n.set)/n.set)
  if(byLC) {
    if(with(grid.sum, n_distinct(p.j.OpI) == n_distinct(p.j.Oth) &
            n_distinct(p.j.Oth) == n_distinct(p.j.Dec) &
            n_distinct(p.j.Dec) == n_distinct(p.j.WP) &
            n_distinct(p.j.WP) == n_distinct(p.j.Evg) &
            n_distinct(p.j.Evg) == n_distinct(p.j.Mxd))) {
      n.LC <- n_distinct(grid.sum$p.j.OpI)
      LC.col <- seq_gradient_pal(low="#e5f5e0", high="#00441b")((1:n.LC)/n.LC)
      LC <- c("OpI", "Oth", "Dec", "WP", "Evg", "Mxd")
    } else {
      cat("DEAL WITH THIS LATER -- need a list for LC.col")
    }
  }
  
  # filenames
  f.nm <- c("pOcc_ad_mn", "pOcc_ad_sd", "pOcc_sb_mn", "pOcc_sb_sd", 
              "pL5_mn", "pL5_sd", "pK_Occ_mn", "pK_Occ_sd")
  f.full <- paste0(p.wd, f.nm)
  
  if(byLC) {
    
    for(l in 1:6) {
      p.mn <- ggplot(grid.sum, aes(x=year, group=p.j)) + ylim(0,100) + 
        labs(x="Year", y="Mean across simulations") +
        scale_colour_manual(paste0(LC[l], ": ", p), values=LC.col) +
        geom_line(aes_string(colour=paste0("p.j.", LC[l])))
      p.sd <- ggplot(grid.sum, aes(x=year, group=p.j)) + 
        labs(x="Year", y="Standard deviation across simulations") +
        scale_colour_manual(paste0(LC[l], ": ", p), values=LC.col) +
        geom_line(aes_string(colour=paste0("p.j.", LC[l])))
      # Adult occupancy
      ggsave(paste0(f.full[1], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.mn + aes_string(y=f.nm[1]) + 
               ggtitle("% occupied cells (adults)"))
      ggsave(paste0(f.full[2], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.sd + aes_string(y=f.nm[2]) + 
               ggtitle("% occupied cells (adults)"))
      # Seed occupancy
      ggsave(paste0(f.full[3], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.mn + aes_string(y=f.nm[3]) + 
               ggtitle("% occupied cells (seeds)"))
      ggsave(paste0(f.full[4], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.sd + aes_string(y=f.nm[4]) + 
               ggtitle("% occupied cells (seeds)"))
      # Low density cells
      ggsave(paste0(f.full[5], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.mn + aes_string(y=f.nm[5]) + 
               ggtitle("% cells with ≤ 5 adults"))
      ggsave(paste0(f.full[6], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.sd + aes_string(y=f.nm[6]) + 
               ggtitle("% cells with ≤ 5 adults"))
      # Occupied cells at K
      ggsave(paste0(f.full[7], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.mn + aes_string(y=f.nm[7]) + 
               ggtitle("% occupied cells at K"))
      ggsave(paste0(f.full[8], "_", LC[l], ".jpg"),  width=w, height=h,
             plot=p.sd + aes_string(y=f.nm[8]) + 
               ggtitle("% occupied cells at K"))
    }
  } else {
    p.mn <- ggplot(grid.sum, aes(x=year, group=p.j)) + ylim(0,100) + 
      labs(x="Year", y="Mean across simulations") +
      scale_colour_manual(p, values=p.col) + 
      geom_line(aes(colour=p.j))
    p.sd <- ggplot(grid.sum, aes(x=year, group=p.j)) + 
      labs(x="Year", y="Standard deviation across simulations") +
      scale_colour_manual(p, values=p.col) + 
      geom_line(aes(colour=p.j))
    
    # Adult occupancy
    ggsave(paste0(f.full[1], ".jpg"), width=w, height=h,
           plot=p.mn + aes_string(y=f.nm[1]) + 
             ggtitle("% occupied cells (adults)"))
    ggsave(paste0(f.full[2], ".jpg"), width=w, height=h,
           plot=p.sd + aes_string(y=f.nm[2]) +
             ggtitle("% occupied cells (adults)"))
    # Seed occupancy
    ggsave(paste0(f.full[3], ".jpg"), width=w, height=h,
           plot=p.mn + aes_string(y=f.nm[3]) +
             ggtitle("% occupied cells (seeds)"))
    ggsave(paste0(f.full[4], ".jpg"), width=w, height=h,
           plot=p.sd + aes_string(y=f.nm[4]) +
             ggtitle("% occupied cells (seeds)"))
    # Low density cells
    ggsave(paste0(f.full[5], ".jpg"), width=w, height=h,
           plot=p.mn + aes_string(y=f.nm[5]) +
             ggtitle("% cells with ≤ 5 adults"))
    ggsave(paste0(f.full[6], ".jpg"), width=w, height=h,
           plot=p.sd + aes_string(y=f.nm[6]) +
             ggtitle("% cells with ≤ 5 adults"))
    # Occupied cells at K
    ggsave(paste0(f.full[7], ".jpg"), width=w, height=h,
           plot=p.mn + aes_string(y=f.nm[7]) +
             ggtitle("% occupied cells at K"))
    ggsave(paste0(f.full[8], ".jpg"), width=w, height=h,
           plot=p.sd + aes_string(y=f.nm[8]) +
             ggtitle("% occupied cells at K"))
  }
  
  # Check for success
  for(f in 1:length(f.nm)) {
    if(sum(grepl(f.nm[f], list.files(p.wd))) > 0) { 
      cat(" ", f.nm[f], "saved\n")
    } else { 
      cat("--- Error:", f.nm[f], "not found! \n") }
  }
}




