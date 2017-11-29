#' Save global and control parameters as RData files
#'
#' This function saves both the global and buckthorn control parameters used in
#' a particular simulation run as a \code{.RData} file in the same directory as
#' the associated plots.
#' @param p.wd Directory to store file in
#' @param age.i Age or age range at first fruiting (e.g., \code{3} or
#'   \code{3-5})
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param control.p NULL or named list of buckthorn control treatment parameters
#'   set with \code{\link{set_control_p}}
#' @return None
#' @keywords parameters, store, save, output
#' @export

save_pars <- function(p.wd, age.i, g.p, control.p) {
  f <- paste0(p.wd, "pars_", age.i, ".RData")
  save(g.p, control.p, file=f)
  if(file.exists(f)) { 
    cat("Parameters stored in", f, "\n") 
  } else {
    cat("--- Error: Parameters not stored \n")
  }
}




#' Save plots of final buckthorn distribution
#'
#' This function saves plots of the final distribution of buckthorn. This
#' includes plots for presence/absence and abundance for both seeds (log
#' abundance) and adults.
#' @param p.wd Directory to store file in
#' @param age.i Age or age range at first fruiting (e.g., \code{3} or
#'   \code{3-5})
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param N.final Dataframe or tibble with processed output from
#'   \code{\link{run_sim}}, including all columns from \code{lc.df} in addition
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.less5}, filtered
#'   to include only the final time step
#' @return None
#' @keywords plots, store, save, output
#' @export

make_plots_final_t <- function(p.wd, age.i, g.p, N.final) {
  require(ggplot2)
  
  # filenames
  ad.ab.f <- paste0(p.wd, "Final_Ab_", age.i, ".jpg")
  sb.ab.f <- paste0(p.wd, "Final_SB_Ab_", age.i, ".jpg")
  ad.pa.f <- paste0(p.wd, "Final_PA_", age.i, ".jpg")
  sb.pa.f <- paste0(p.wd, "Final_SB_PA_", age.i, ".jpg")
  
  # Adult abundance
  ad.ab.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.adult, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste("Adult abundance. Year", g.p$tmax+1))
  ggsave(ad.ab.f, ad.ab.fin, width=8, height=6, units="in")
  rm(ad.ab.fin)
  
  # Seed bank log abundance
  sb.ab.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.sb, colour=inbd)) +
    geom_tile() + scale_fill_gradient(low="white", high="red") +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste("log seed abundance. Year", g.p$tmax+1))
  ggsave(sb.ab.f, sb.ab.fin, width=8, height=6, units="in")
  rm(sb.ab.fin)
  
  # Adult presence/absence
  adult.pa.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.adult>0, colour=inbd)) +
    geom_tile() + scale_fill_manual(values=c("FALSE"="white", "TRUE"="red")) + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste("Adult presence. Year", g.p$tmax+1))
  ggsave(ad.pa.f, adult.pa.fin, width=8, height=6, units="in")
  rm(adult.pa.fin)
  
  # Seed bank presence/absence
  sb.pa.fin <- ggplot(N.final, aes(x=x, y=-y, fill=N.sb>0, colour=inbd)) +
    geom_tile() + scale_fill_manual(values=c("FALSE"="white", "TRUE"="red")) +
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle(paste("Seed presence. Year", g.p$tmax+1))
  ggsave(sb.pa.f, sb.pa.fin, width=8, height=6, units="in")
  rm(sb.pa.fin)
  
  # Check for success
  if(file.exists(ad.ab.f)) { 
    cat("Final adult abundance saved as", ad.ab.f, "\n") 
  } else {
    cat("--- Error: Final adult abundance not saved! \n")
  }
  if(file.exists(sb.ab.f)) { 
    cat("Final seed bank abundance saved as", sb.ab.f, "\n") 
  } else {
    cat("--- Error: Final seed bank abundance not saved! \n")
  }
  if(file.exists(ad.pa.f)) { 
    cat("Final adult presence saved as", ad.pa.f, "\n") 
  } else {
    cat("--- Error: Final adult presence not saved! \n")
  }
  if(file.exists(sb.pa.f)) { 
    cat("Final seed bank presence saved as", sb.pa.f, "\n") 
  } else {
    cat("--- Error: Final seed bank presence not saved! \n")
  }
}




#' Save gifs of buckthorn distribution through time
#'
#' This function saves gifs of the distribution of buckthorn. This includes
#' plots for presence/absence and abundance for both seeds (log abundance) and
#' adults. There is also a gif distinguishing cells with a low density of
#' adults.
#' @param p.wd Directory to store file in
#' @param age.i Age or age range at first fruiting (e.g., \code{3} or
#'   \code{3-5})
#' @param g.p Named list of global parameters set with \code{\link{set_g_p}}
#' @param N.out Dataframe or tibble with processed output from
#'   \code{\link{run_sim}}, including all columns from \code{lc.df} in addition
#'   to \code{year}, \code{N.adult}, \code{N.sb}, and \code{N.less5}
#' @return None
#' @keywords plots, gif, store, save, output
#' @export

make_plots_gifs <- function(p.wd, age.i, g.p, N.out) {
  require(gganimate)
  
  # filenames
  ad.ab.f <- paste0(p.wd, "Ab_", age.i, ".gif")
  sb.ab.f <- paste0(p.wd, "SB_Ab_", age.i, ".gif")
  ad.pa.f <- paste0(p.wd, "PA_", age.i, ".gif")
  sb.pa.f <- paste0(p.wd, "SB_PA_", age.i, ".gif")
  ad.lo.f <- paste0(p.wd, "LoDens_", age.i, ".gif")
  
  # Adult abundance
  ad.ab <- ggplot(N.out, aes(x=x, y=-y, fill=N.adult, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle("Adult abundance. Year")
  gganimate(ad.ab, ad.ab.f, interval=0.2, ani.width=800, ani.height=600)
  rm(ad.ab)
  
  # Seed bank log abundance
  sb.ab <- ggplot(N.out, aes(x=x, y=-y, fill=N.sb, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_gradient(low="white", high="red") + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle("log seed abundance. Year")
  gganimate(sb.ab, sb.ab.f, interval=0.2, ani.width=800, ani.height=600)
  rm(sb.ab)
  
  # Adult presence/absence
  adult.pa <- ggplot(N.out, aes(x=x, y=-y, fill=N.adult>0, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_manual(values=c("white","red")) + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle("Adult presence. Year")
  gganimate(adult.pa, ad.pa.f, interval=0.2, ani.width=800, ani.height=600)
  rm(adult.pa)
  
  # Seed bank presence/absence
  sb.pa <- ggplot(N.out, aes(x=x, y=-y, fill=N.sb > 0, frame=year, colour=inbd)) + 
    geom_tile() + scale_fill_manual(values=c("white","red")) + 
    scale_colour_manual(values=c("gray", NA)) + 
    ggtitle("Seed presence. Year")
  gganimate(sb.pa, sb.pa.f, interval=0.2, ani.width=800, ani.height=600)
  rm(sb.pa)
  
  # Low adult density
  adult.lo <- ggplot(N.out, aes(x=x, y=-y, fill=N.less5, frame=year)) +
    geom_tile() + ggtitle("Blue: Adult density ≤ 5. Year") + 
    scale_fill_gradient2(low="white", mid="blue", high="pink", midpoint=1)
  gganimate(adult.lo, ad.lo.f, interval=0.2, ani.width=800, ani.height=600)
  rm(adult.lo)
  
  # Check for success
  if(file.exists(ad.ab.f)) { 
    cat("Adult abundance saved as", ad.ab.f, "\n") 
  } else {
    cat("--- Error: Adult abundance not saved! \n")
  }
  if(file.exists(sb.ab.f)) { 
    cat("Seed bank abundance saved as", sb.ab.f, "\n") 
  } else {
    cat("--- Error: Seed bank abundance not saved! \n")
  }
  if(file.exists(ad.pa.f)) { 
    cat("Adult presence saved as", ad.pa.f, "\n") 
  } else {
    cat("--- Error: Adult presence not saved! \n")
  }
  if(file.exists(sb.pa.f)) { 
    cat("Seed bank presence saved as", sb.pa.f, "\n") 
  } else {
    cat("--- Error: Seed bank presence not saved! \n")
  }
  if(file.exists(ad.lo.f)) { 
    cat("Low adult abundance saved as", ad.lo.f, "\n") 
  } else {
    cat("--- Error: Low adult abundance not saved! \n")
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
  require(ggplot2)
  
  # filenames
  opi.f <- paste0(p.wd, "LC_OpI.jpg")
  oth.f <- paste0(p.wd, "LC_Oth.jpg")
  dec.f <- paste0(p.wd, "LC_Dec.jpg")
  wp.f <- paste0(p.wd, "LC_WP.jpg")
  evg.f <- paste0(p.wd, "LC_Evg.jpg")
  mxd.f <- paste0(p.wd, "LC_Mxd.jpg")
  
  p.opi <- ggplot(lc.df, aes(x=x, y=-y, fill=OpI, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="red", limits=c(0,1))
  ggsave(opi.f, p.opi, width=10, height=7)
  rm(p.opi)
  
  p.oth <- ggplot(lc.df, aes(x=x, y=-y, fill=Oth, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="gray30", limits=c(0,1))
  ggsave(oth.f, p.oth, width=10, height=7)
  rm(p.oth)
  
  p.dec <- ggplot(lc.df, aes(x=x, y=-y, fill=Dec, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="green3", limits=c(0,1))
  ggsave(dec.f, p.dec, width=10, height=7)
  rm(p.dec)
  
  p.wp <- ggplot(lc.df, aes(x=x, y=-y, fill=WP, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="orchid", limits=c(0,1))
  ggsave(wp.f, p.wp, width=10, height=7)
  rm(p.wp)
  
  p.evg <- ggplot(lc.df, aes(x=x, y=-y, fill=Evg, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="darkgreen", limits=c(0,1))
  ggsave(evg.f, p.evg, width=10, height=7)
  rm(p.evg)
  
  p.mxd <- ggplot(lc.df, aes(x=x, y=-y, fill=Mxd, colour=inbd)) + geom_tile() +
    scale_colour_manual(values=c("gray", NA)) + 
    scale_fill_gradient(low="white", high="yellowgreen", limits=c(0,1))
  ggsave(mxd.f, p.mxd, width=10, height=7)
  rm(p.mxd)
  
  # Check for success
  if(file.exists(opi.f)) { 
    cat("Open invasible proportions saved as", opi.f, "\n") 
  } else {
    cat("--- Error: Open invasible proportions not saved! \n")
  }
  if(file.exists(oth.f)) { 
    cat("Other proportions saved as", oth.f, "\n") 
  } else {
    cat("--- Error: Other proportions not saved! \n")
  }
  if(file.exists(dec.f)) { 
    cat("Deciduous forest proportions saved as", dec.f, "\n") 
  } else {
    cat("--- Error: Deciduous forest proportions not saved! \n")
  }
  if(file.exists(wp.f)) { 
    cat("White pine forest proportions saved as", wp.f, "\n") 
  } else {
    cat("--- Error: White pine forest proportions not saved! \n")
  }
  if(file.exists(evg.f)) { 
    cat("Evergreen forest proportions saved as", evg.f, "\n") 
  } else {
    cat("--- Error: Evergreen forest proportions not saved! \n")
  }
  if(file.exists(mxd.f)) { 
    cat("Mixed forest proportions saved as", mxd.f, "\n") 
  } else {
    cat("--- Error: Mixed forest proportions not saved! \n")
  }
}






