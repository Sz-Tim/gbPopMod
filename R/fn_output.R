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
  library(ggplot2); options(bitmapType='cairo')
  theme_set(theme_bw() + theme(panel.grid=element_blank()))
  
  # filenames
  f.nm <- c("ad_Ab", "sb_Ab", "ad_sd", "sb_sd", "ad_pP", "sb_pP")
  f.full <- paste0(p.wd, "Final_", f.nm, ".jpg") 
  
  p.fin <- ggplot(N.final, aes(x=lon, y=lat)) +
    scale_fill_gradient(low="white", high="red") +
    theme(panel.background=element_rect(fill="gray30"))
  
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
           ggtitle(paste(txt, "Adult pr(presence). Year", g.p$tmax+1)))
  ggsave(f.full[6], width=w, height=h, 
         plot=p.fin + geom_tile(aes_string(fill=f.nm[6])) +
           ggtitle(paste(txt, "Seed pr(presence). Year", g.p$tmax+1)))
  
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
#' @return None
#' @keywords plots, gif, store, save, output
#' @export

make_plots_gifs <- function(p.wd, g.p, N.out, txt=NULL, w=800, h=600) {
  library(gganimate); options(bitmapType='cairo')
  theme_set(theme_bw() + theme(panel.grid=element_blank()))
  
  # filenames
  f.nm <- c("ad_Ab", "sb_Ab", "ad_sd", "sb_sd", "ad_pP", "sb_pP", "ad_L5")
  f.full <- paste0(p.wd, f.nm, ".gif")
  
  p.gif <- ggplot(N.out, aes(x=lon, y=lat)) + 
    scale_fill_gradient(low="white", high="red") + 
    theme(panel.background=element_rect(fill="gray30")) +
    transition_time(as.numeric(year))
  
  # gifs
  anim_save(f.full[1],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[1])) + 
              ggtitle(paste(txt, "Adult mean(abundance). Year {frame_time}")),
              nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[2],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[2])) + 
                ggtitle(paste(txt, "Seed mean(log abundance). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[3],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[3])) + 
                ggtitle(paste(txt, "Adult sd(abundance). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[4],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[4])) + 
                ggtitle(paste(txt, "Seed sd(log abundance). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[5],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[5])) + 
                ggtitle(paste(txt, "Adult pr(presence). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[6],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[6])) + 
                ggtitle(paste(txt, "Seed pr(presence). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  anim_save(f.full[7],
            animate(p.gif + geom_tile(aes_string(fill=f.nm[7])) + 
                scale_fill_gradient2(low="white", mid="blue", 
                                      high="pink", midpoint=1) + 
                ggtitle(paste(txt, "Adult mean(0 < N < 5). Year {frame_time}")),
                nframes=n_distinct(N.out$year), width=w, height=h, units="px"))
  
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
  library(ggplot2); options(bitmapType='cairo')
  theme_set(theme_bw() + theme(panel.grid=element_blank()))
  
  # filenames
  f.full <- paste0(p.wd, c("LC_Opn.jpg", "LC_Oth.jpg", "LC_Dec.jpg", 
                           "LC_WP.jpg", "LC_Evg.jpg", "LC_Mxd.jpg"))
  
  p.lc <- ggplot(lc.df, aes(x=lon, y=lat, colour=inbd)) + 
    theme(panel.background=element_rect(fill="gray30")) +
    scale_colour_manual(values=c("gray", NA))
  
  ggsave(f.full[1], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Opn)) + ggtitle("Open Invasible") +
           scale_fill_gradient(low="white", high="red", limits=c(0,1)))
  ggsave(f.full[2], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Oth)) + ggtitle("Other")+
           scale_fill_gradient(low="white", high="gray30", limits=c(0,1)))
  ggsave(f.full[3], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Dec)) + ggtitle("Deciduous Forest")+
           scale_fill_gradient(low="white", high="green3", limits=c(0,1)))
  ggsave(f.full[4], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=WP)) + ggtitle("White Pine Forest")+
           scale_fill_gradient(low="white", high="orchid", limits=c(0,1)))
  ggsave(f.full[5], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Evg)) + ggtitle("Evergreen Forest")+
           scale_fill_gradient(low="white", high="darkgreen", limits=c(0,1)))
  ggsave(f.full[6], width=w, height=h,
         plot=p.lc + geom_tile(aes(fill=Mxd)) + ggtitle("Mixed Forest")+
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
#' @param cell.sum Dataframe or tibble with processed output for summaries of
#'   time-to-event by cell
#' @param byLC \code{FALSE} For LC-specific parameters, should separate plots be
#'   created for each LC?
#' @param txt Text to append to start of plot title
#' @param w Figure output width (inches)
#' @param h Figure output height (inches)
#' @return None
#' @keywords plots, store, save, output
#' @export

make_plots_gridSummary <- function(p, p.wd, grid.sum, cell.sum, byLC=FALSE, 
                                   txt=NULL, w=8, h=6) {
  library(scales); library(ggplot2)
  theme_set(theme_bw() + theme(panel.grid=element_blank()))
  options(bitmapType='cairo')
  
  n.set <- length(unique(grid.sum$p.j))
  p.col <- seq_gradient_pal(low="#e5f5e0", high="#00441b")((1:n.set)/n.set)
  
  # filenames
  f.nm <- c("pOcc_ad_mn", "pOcc_ad_sd", "pOcc_sb_mn", "pOcc_sb_sd", 
            "pL5_mn", "pL5_sd", "pK_mn", "pK_sd", "pK_Occ_mn", "pK_Occ_sd",
            "t0K_mn", "t0K_sd", "tL5_mn", "tL5_sd")
  f.full <- paste0(p.wd, f.nm)
  p.mn <- ggplot(grid.sum, aes(x=year, group=p.j, colour=p.j)) + 
    labs(x="Year") + ylim(0,100) + 
    scale_colour_manual(p[1], values=p.col) + geom_line()
  p.sd <- ggplot(grid.sum, aes(x=year, group=p.j, colour=p.j)) + 
    labs(x="Year") +
    scale_colour_manual(p[1], values=p.col) + geom_line()
  t.mn <- ggplot(cell.sum, aes(group=p.j, colour=p.j)) +  
    labs(y="Density") +
    scale_colour_manual(p[1], values=p.col) + geom_density()
  t.sd <- ggplot(cell.sum, aes(group=p.j, colour=p.j)) +  
    labs(y="Density") +
    scale_colour_manual(p[1], values=p.col) + geom_density()
  
  # Adult occupancy
  ggsave(paste0(f.full[1], ".jpg"), width=w, height=h,
         plot=p.mn + aes_string(y=f.nm[1]) + 
           labs(title="% occupied cells (adults)", 
                y="Mean across simulations"))
  ggsave(paste0(f.full[2], ".jpg"), width=w, height=h,
         plot=p.sd + aes_string(y=f.nm[2]) +
           labs(title="% occupied cells (adults)", 
                y="Std dev across simulations"))
  # Seed occupancy
  ggsave(paste0(f.full[3], ".jpg"), width=w, height=h,
         plot=p.mn + aes_string(y=f.nm[3]) +
           labs(title="% occupied cells (seeds)",
                y="Mean across simulations"))
  ggsave(paste0(f.full[4], ".jpg"), width=w, height=h,
         plot=p.sd + aes_string(y=f.nm[4]) +
           labs(title="% occupied cells (seeds)", 
                y="Std dev across simulations"))
  # Low density cells
  ggsave(paste0(f.full[5], ".jpg"), width=w, height=h,
         plot=p.mn + aes_string(y=f.nm[5]) +
           labs(title="% cells with ≤ 5 adults",
                y="Mean across simulations"))
  ggsave(paste0(f.full[6], ".jpg"), width=w, height=h,
         plot=p.sd + aes_string(y=f.nm[6]) +
           labs(title="% cells with ≤ 5 adults", 
                y="Std dev across simulations"))
  # Cells at K
  ggsave(paste0(f.full[7], ".jpg"), width=w, height=h,
         plot=p.mn + aes_string(y=f.nm[7]) +
           labs(title="% cells at K", 
                y="Mean across simulations"))
  ggsave(paste0(f.full[8], ".jpg"), width=w, height=h,
         plot=p.sd + aes_string(y=f.nm[8]) +
           labs(title="% cells at K", 
                y="Std dev across simulations"))
  # Occupied cells at K
  ggsave(paste0(f.full[9], ".jpg"), width=w, height=h,
         plot=p.mn + aes_string(y=f.nm[9]) +
           labs(title="% occupied cells at K", 
                y="Mean across simulations"))
  ggsave(paste0(f.full[10], ".jpg"), width=w, height=h,
         plot=p.sd + aes_string(y=f.nm[10]) +
           labs(title="% occupied cells at K", 
                y="Std dev across simulations"))
  # Time from N=1 to N=K
  ggsave(paste0(f.full[11], ".jpg"), width=w, height=h,
         plot=t.mn + aes_string(x=f.nm[11]) + 
           labs(title="Time from N=1 to N=K", x="mean(years)"))
  ggsave(paste0(f.full[12], ".jpg"), width=w, height=h,
         plot=t.sd + aes_string(x=f.nm[12]) + 
           labs(title="Time from N=1 to N=K", x="sd(years)"))
  # Time from N=1 to N>5
  ggsave(paste0(f.full[13], ".jpg"), width=w, height=h,
         plot=t.mn + aes_string(x=f.nm[13]) + 
           labs(title="Time from N=1 to N>5", x="mean(years)"))
  ggsave(paste0(f.full[14], ".jpg"), width=w, height=h,
         plot=t.sd + aes_string(x=f.nm[14]) + 
           labs(title="Time from N=1 to N>5", x="sd(years)"))    
  
  # Check for success
  for(f in 1:length(f.nm)) {
    if(sum(grepl(f.nm[f], list.files(p.wd))) > 0) { 
      cat(" ", f.nm[f], "saved\n")
    } else { 
      cat("--- Error:", f.nm[f], "not found! \n") }
  }
}




