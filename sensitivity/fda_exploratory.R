# Functional Data Analysis Exploration

# This script explores the possibility for functional data analysis as a method
# of interpreting the output from the sensitivity analyses. Because of the
# model structure, the population will ultimately converge to K in each cell, 
# but the parameter changes will affect how and how quickly the final state is
# reached. Thus, we are most interested in the shapes of the curves rather than
# any single summary statistic. 


# set environment
Packages <- c("gbPopMod", "tidyverse", "magrittr", "stringr", "here", "doSNOW",
              "fastmatch", "scales", "gganimate", "fda")
suppressMessages(invisible(lapply(Packages, library, character.only=TRUE)))
run_wrap <- TRUE
if(run_wrap) source("hpc/hpc_wrap.R")


# focus on one parameter
p <- readRDS("hpc/p.rds")[6]
p.seq <- list(readRDS("hpc/p_seq.rds")[6])
par.wd <- paste0("out/327_t50/", p)
parSet.wd <- list.dirs(par.wd, recursive=FALSE)

grid.sum <- map_df(parSet.wd, 
                   ~suppressMessages(read_csv(paste0(., "/grid_j.csv")))) %>%
  mutate(year=as.numeric(year)) 
cell.sum <- map_df(parSet.wd, 
                   ~suppressMessages(read_csv(paste0(., "/cell_j.csv"))))
if(length(p.seq[[1]])==1) {
  grid.sum %<>% mutate(p.j=as.factor(p.j)) 
  cell.sum %<>% mutate(p.j=as.factor(p.j)) 
} else {
  grid.sum %<>% mutate_at(vars(p.j:p.j.Mxd), as.factor)
  cell.sum %<>% mutate_at(vars(p.j:p.j.Mxd), as.factor)
}

grd <- list(year=grid.sum$year[1:51], 
            pOcc_sb_mn=matrix(grid.sum$pOcc_sb_mn, nrow=51))
grd.basis <- with(grd, smooth.basisPar(argvals=year, y=pOcc_sb_mn, lambda=.001))
plot(grd.basis$fd, xlab="Year", ylab="Mean Occupancy", ylim=c(0,100), 
     lty=1, col=topo.colors(length(parSet.wd)+5)[1:length(parSet.wd)])
plot(deriv(grd.basis$fd), xlab="Year", ylab="Rate of spread (% cells/yr)",
     lty=1, col=topo.colors(length(parSet.wd)+5)[1:length(parSet.wd)])


t0K_mn.mx <- list(x=NULL, y=NULL)
for(j in 1:length(parSet.wd)) {
  j.dens <- density(filter(cell.sum, p.j==p.seq[[1]][[1]][j])$t0K_mn,
                    from=0, to=50)
  t0K_mn.mx$x <- j.dens$x
  t0K_mn.mx$y <- cbind(t0K_mn.mx$y, j.dens$y)
}

cll.basis <- with(t0K_mn.mx, smooth.basisPar(argvals=x, y=y, lambda=0.1))
plot(cll.basis$fd, xlab="mean(years)", ylab="density", 
     lty=1, col=topo.colors(length(parSet.wd)+5)[1:length(parSet.wd)])
plot(deriv(cll.basis$fd), xlab="mn(years)", ylab="Rate of spread (% cells/yr)",
     lty=1, col=topo.colors(length(parSet.wd)+5)[1:length(parSet.wd)])




