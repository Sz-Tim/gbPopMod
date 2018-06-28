#' Antilogit function
#'
#' Calculate the antilogit for the object \code{x}, returning the corresponding
#' probability (range 0-1)
#' @param x Vector of unconstrained values
#' @return Vector of constrained values
#' @keywords antilogit
#' @export

antilogit <- function (x) {
  exp(x)/(1 + exp(x))
}




#' Logit function
#'
#' Calculate the logit for the probability vector \code{x}, returning the
#' corresponding unconstrained value
#' @param x Vector of constrained values
#' @return Vector of unconstrained values
#' @keywords antilogit
#' @export

logit <- function (x) {
  log(x/(1-x))
}




#' Expand all pairwise combinations of two vectors into one character vector
#'
#' Similar to \code{\link[base]{expand.grid}} but inputs two vectors and returns
#' a single character vector with the values from the two vectors separated by
#' "_" by default.
#' @param x Vector
#' @param y Vector
#' @param sep Like paste, specifies character(s) to separate vector values
#' @return Vector (character)
#' @keywords expand.grid
#' @export

expand_v <- function(x, y, sep="_") {
  paste(rep.int(x, length(y)), 
        rep.int(y, rep.int(length(x),length(y))),
        sep=sep)
}




#' Expand land cover-based parameter ranges
#'
#' Similar to \code{\link[base]{expand.grid}} but inputs a vector of minimum
#' values and a vector of maximum values in addition to a length.out parameter.
#' It returns a list of vectors, with an element for each combination of land
#' cover parameters. By default, \code{all.combo=FALSE} and all land covers are
#' incremented jointly. If \code{all.combo=TRUE}, then all combinations of land
#' cover parameter values will be output. If \code{LC="all"}, then all land
#' cover types will be incremented across the specifie range. If \code{LC} is
#' set to a specific land cover category, then only that category will be
#' incremented while the other categories are held constant.
#' @param gp Named list of global parameters. If \code{LC != "all"}, then
#'   default values from g.p are used for the land cover categories that are not
#'   being varied.
#' @param param \code{NULL} Character scalar of which parameter to vary. If
#'   \code{LC != "all"}, then this must be specified.
#' @param LC Character scalar of which land cover categories to increment; must
#'   take one of \code{c("all", "Opn", "Oth", "Dec", "WP", "Evg", "Mxd")}
#' @param all.combo \code{FALSE} Should all combinations of all land cover
#'   values be generated?
#' @param len.out Number of parameter values per land cover category
#' @param lc.min Vector \code{length=n.lc} of minimum parameter values.
#' @param lc.max Vector \code{length=n.lc} of maximum parameter values
#' @return List of vectors, each length 6
#' @keywords expand.grid, sensitivity
#' @export

expand_LCs <- function(gp=g.p, param=NULL, LC="all", all.combo=FALSE, 
                       len_out=6, lc.min=rep(0.1, 6), lc.max=rep(0.9, 6)) {
  library(tidyverse); library(purrr)
  names(lc.min) <- c("Opn", "Oth", "Dec", "WP", "Evg", "Mxd")
  names(lc.max) <- names(lc.min)
  if(LC=="all") {
    g <- map2_df(lc.min, lc.max, ~seq(.x, .y, length.out=len_out))
    if(all.combo) {
      g <- g %>% expand.grid %>% as.matrix
    } else {
      g <- as.matrix(g)
    }
  } else {
    g <- t(replicate(len_out, gp[[param]]))
    colnames(g) <- names(lc.min)
    g[,LC] <- seq(lc.min[LC], lc.max[LC], length.out=len_out)
  } 
  return(l=lapply(seq_len(nrow(g)), function(i) g[i,]))
}




#' Expand to all combinations of canopy parameter ranges
#'
#' Similar to \code{\link[base]{expand.grid}} but inputs two vectors with min
#' and max parameter values (one set for open canopy and one set for closed
#' canopy) in addition to a length.out parameter, and returns a list of vectors,
#' with an element for each combination of land cover parameters.
#' @param Open Vector of min & max for open canopy categories (open invasible
#'   and other)
#' @param Closed Vector of min & max for closed canopy categories (deciduous,
#'   white pine, other evergreen, and mixed forests)
#' @param length.out Number of parameter values per land cover category
#' @return List of vectors, each length 6
#' @keywords expand.grid, sensitivity
#' @export

expand_cnpy <- function(Op=c(3, 7), Cl=c(3, 7), length_out=2) {
  g <- expand.grid(Op=seq(Op[1], Op[2], length.out=length_out),
                   Cl=seq(Cl[1], Cl[2], length.out=length_out))
  g <- t(apply(g, 1, function(i) rep(i, times=c(2,4))))
  return(l=lapply(seq_len(nrow(g)), function(i) g[i,]))
}




#' Aggregate compositional data within each cell
#'
#' Reformat and calculate expected cell-means based on land cover composition
#' for relevant parameters.
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions,
#'   other covariates, and cell id info
#' @param K Vector \code{length=n.lc} with carrying capacity for each land cover
#'   type or vector of slopes corresponding with columns in lc.df
#' @param s.M Vector \code{length=n.lc} with juvenile survival probability for
#'   each land cover type or vector of slopes corresponding with columns in
#'   lc.df
#' @param s.N \code{c(1, 1, 1, 1, 1, 1)} Vector \code{length=n.lc} of annual
#'   adult survival rates or vector of slopes corresponding with columns in
#'   lc.df
#' @param mu Vector \code{length=n.lc} with mean per-individual fruit production
#'   for each land cover type or vector of slopes corresponding with columns in
#'   lc.df
#' @param p.f Vector \code{length=n.lc} with mean probability of fruiting for
#'   each land cover type or vector of slopes corresponding with columns in
#'   lc.df
#' @param p.c Vector \code{length=n.lc} with proportion of fruits eaten by birds
#'   for each land cover type or vector of slopes corresponding with columns in
#'   lc.df
#' @param p Vector \code{length=n.lc} with seedling establishment probability
#'   for each land cover type or vector of slopes corresponding with columns in
#'   lc.df
#' @param p.trt Tibble with grid id and modified establishment probabilities for
#'   cells with ground cover treatments; default = NULL
#' @param edges Character taking the value of one of: \code{"wall", "sink",
#'   "none"} where \code{"wall"} results in a dispersal probability of 0 for all
#'   out-of-bound cells with no populations modeled, \code{"sink"} results in
#'   dispersal of seeds to out-of-bound cells but no populations modeled, and
#'   \code{"none"} results in dispersal of seeds and populations modeled
#' @param method \code{"wt.mn"} Method for calculating cell expectations, taking
#'   values of \code{"wt.mn"} or \code{"lm"}. If \code{"wt.mn"}, the expectation
#'   for each parameter is the weighted mean across land cover types
#'   proportional to their coverage, with the land cover specific values stored
#'   in the parameter vectors. If \code{"lm"}, the expectation is calculated in
#'   a regression with the slopes contained in each parameter vector.
#'   Individuals cannot be assigned to specific land cover categories with
#'   \code{"lm"}, so \code{"m"} must be scalar.
#' @return Named list with values aggregated within cells based on land cover
#'   types. Includes: \describe{ \item{\code{lc.mx}}{Matrix \code{(ncol=n.lc,
#'   nrow=ngrid)} with land cover proportions} \item{\code{K.E}}{Vector
#'   \code{length=ngrid} with total K} \item{\code{K.lc}}{Matrix
#'   \code{(ncol=n.lc, nrow=ngrid)} with K per land cover category}
#'   \item{\code{s.M.E}}{Vector \code{length=ngrid} with pr(juvenile surv)}
#'   \item{\code{rel.dens}}{Matrix \code{(ncol=n.lc, nrow=ngrid)} with relative
#'   density among land cover categories} \item{\code{mu.E}}{Vector
#'   \code{length=ngrid} with mean fruit produced per adult)}
#'   \item{\code{p.f.E}}{Vector \code{length=ngrid} with fruiting probability}
#'   \item{\code{p.c.E}}{Vector \code{length=ngrid} with proportion eaten by
#'   birds} \item{\code{p.E}}{Vector \code{length=ngrid} with seedling
#'   establishment probabilities} }
#' @note If \code{method="lm"}, then each parameter vector will be treated as a
#'   set of slopes for the covariates in lc.df with the number of covariates
#'   used in each regression is \code{length(param)-1} and the first element of
#'   \code{param} is the intercept.
#' @note If \code{!is.null(p.trt)}, then the associated p.E values are
#'   substituted in the cells that received a relevant management treatments.
#' @keywords premultiply, aggregate, set up, initialize
#' @export

cell_E <- function(lc.df, K, s.M, s.N, mu, p.f, p.c, p, 
                   p.trt=NULL, edges="wall", method="wt.mn") {
  
  if(method=="wt.mn") {
    # scalar = same value for all LC categories
    if(length(K)==1) K <- rep(K, 6)
    if(length(s.N)==1) s.N <- rep(s.N, 6)
    if(length(s.M)==1) s.M <- rep(s.M, 6)
    if(length(mu)==1) mu <- rep(mu, 6)
    if(length(p.f)==1) p.f <- rep(p.f, 6)
    if(length(p.c)==1) p.c <- rep(p.c, 6)
    if(length(p)==1) p <- rep(p, 6)
    # take weighted mean
    lc.mx <- as.matrix(lc.df[,4:9])
    K.E <- round(lc.mx %*% K)
    K.lc <- round(t(t(lc.mx) * K))
    rel.dens <- t(apply(lc.mx, 1, function(x) K*x/c(x%*%K)))
    s.N.E <- c(lc.mx %*% s.N)
    s.M.E <- c(lc.mx %*% s.M)
    mu.E <- lc.mx %*% mu
    p.f.E <- lc.mx %*% p.f
    p.c.E <- lc.mx %*% p.c
    p.E <- lc.mx %*% p
  } else if(method=="lm") {
    lc.mx <- cbind(1, as.matrix(select(lc.df, 
                              -one_of("x", "y", "x_y", "inbd", "id", "id.in"))))
    K.E <- exp(lc.mx[,1:length(K)] %*% K)
    K.lc <- NULL
    rel.dens <- NULL
    s.N.E <- c(antilogit(lc.mx[,1:length(s.N)] %*% s.N))
    s.M.E <- c(antilogit(lc.mx[,1:length(s.M)] %*% s.M))
    mu.E <- exp(lc.mx[,1:length(mu)] %*% mu)
    p.f.E <- antilogit(lc.mx[,1:length(p.f)] %*% p.f)
    p.c.E <- antilogit(lc.mx[,1:length(p.c)] %*% p.c)
    p.E <- antilogit(lc.mx[,1:length(p)] %*% p)
  }
  
  if(!is.null(p.trt)) {
    p.E[p.trt$id,] <- p.trt$p
  }
  
  if(edges=="sink") p.E[!lc.df$inbd] <- 0
  
  return(list(lc.mx=lc.mx, K.E=K.E, K.lc=K.lc, rel.dens=rel.dens,
              s.M.E=s.M.E, s.N.E=s.N.E, mu.E=mu.E, p.f.E=p.f.E,
              p.c.E=p.c.E, p.E=p.E))
}




#' Initialize populations randomly
#'
#' Initialize populations randomly with populated cells containing adults at 50%
#' K and juveniles at 10% K
#' @param ngrid Number of grid cells in entire map
#' @param g.p Named list of global parameters
#' @param lc.df Dataframe or tibble with xy coords, land cover proportions, and
#'   cell id info
#' @return Matrix or array of initial abundances with \code{dim=c(ngrid, (n.lc),
#'   m.max)}
#' @keywords initialize, set up
#' @export

pop_init <- function(ngrid, g.p, lc.df) {
  
  p.0 <- sample(lc.df$id[lc.df$inbd], g.p$N.p.t0)
  m.max <- max(g.p$m)  # adult age bin
  
  if(length(unique(g.p$m)) == 1) {
    N.init <- matrix(0, ngrid, m.max)  # column for each age class
    N.init[p.0,m.max] <- round(as.matrix(lc.df[lc.df$id %in% p.0, 4:9]) %*% 
                                (g.p$K/2))
    N.init[p.0,-m.max] <- round(N.init[p.0,m.max]/5)
    
  } else {
    N.init <- array(0, dim=c(ngrid, g.p$n.lc, m.max))
    N.init[p.0,,m.max] <- round(t(t(as.matrix(lc.df[lc.df$id %in% p.0, 4:9])) * 
                                   g.p$K/2))
    N.init[p.0,,-m.max] <- round(N.init[p.0,,m.max]/5)
  }
  
  return(N.init)
}




#' Set global parameters
#'
#' Set the global parameters for the simulation, allowing for individual
#' elements to be reassigned.
#' @param tmax \code{100} Number of time steps per simulation
#' @param dem.st \code{FALSE} Include stochasticity in demography?
#' @param sdd.st \code{TRUE} Include stochasticity in short distance dispersal?
#' @param bank \code{TRUE} Include seedbank?
#' @param n.cores \code{4} Number of cores for parallelizing sdd.pr calculation
#' @param lc.r \code{100} Maximum number of rows (\code{y}) in landscape
#' @param lc.c \code{100} Maximum number of columns (\code{x}) in landscape
#' @param n.lc \code{6} Number of land cover categories
#' @param N.p.t0 \code{10} Number of cells with buckthorn at t=1
#' @param p.f \code{c(0.9, 0.1, 0.29, 0.23, 0.2, 0.3)} Vector \code{length=n.lc}
#'   of fruiting probabilities
#' @param mu \code{c(200, 100, 40, 20, 20, 10)} Vector \code{length=n.lc} of
#'   mean fruit per adult
#' @param gamma \code{2.3} Scalar: mean number of seeds per fruit
#' @param m \code{4} Vector \code{length=n.lc} or scalar of age at first
#'   fruiting. Individuals at this age are considered adults
#' @param p.c \code{c(0.3, 0.1, 0.2, 0.2, 0.2, 0.1)} Vector \code{length=n.lc}
#'   of proportion of fruits eaten by birds, with \code{1-p.c} assumed to drop
#'   directly below buckthorn individuals
#' @param sdd.rate \code{0.1} 1/mn for exponential dispersal kernel
#' @param sdd.max \code{15} Maximum dispersal distance in cells
#' @param bird.hab \code{c(0.35, 0.35, 0.05, 0.1, 0.1, 0.05)} Vector
#'   \code{length=n.lc} of bird habitat preferences
#' @param n.ldd \code{1} Number of long distance dispersal events per year
#' @param s.c \code{0.6} Seed viability post-digestion
#' @param s.B \code{0.3} Probability of annual survival in seed bank
#' @param s.M \code{c(0.9, 0.1. 0.6, 0.6, 0.6, 0.6)} Vector \code{length=n.lc}
#'   of annual juvenile survival rates
#' @param s.N \code{c(1, 1, 1, 1, 1, 1)} Vector \code{length=n.lc} of annual
#'   adult survival rates
#' @param K \code{c(750, 10, 100, 100, 300, 100)} Vector (length=n.lc) of
#'   carrying capacities for adults
#' @param g.D \code{0.5} Probability of direct germination (i.e., a seed
#'   germinates in the same year it is produced)
#' @param g.B \code{0.5} Probability of germinating from the seed bank
#' @param p \code{c(0.07, 0.01, 0.08, 0.02, 0.02, 0.03)} Vector
#'   \code{length=n.lc} of seedling establishment probabilities
#' @param edges \code{"wall"} Boundary behavior, taking values of \code{"wall"},
#'   \code{"sink"}, or \code{"none"}. See boundary_behavior.Rmd for descriptions
#' @param method \code{"wt.mn"} Method for calculating cell expectations, taking
#'   values of \code{"wt.mn"} or \code{"lm"}. If \code{"wt.mn"}, the expectation
#'   for each parameter is the weighted mean across land cover types
#'   proportional to their coverage, with the land cover specific values stored
#'   in the parameter vectors. If \code{"lm"}, the expectation is calculated in
#'   a regression with the slopes contained in each parameter vector.
#'   Individuals cannot be assigned to specific land cover categories with
#'   \code{"lm"}, so \code{"m"} must be scalar.
#' @return Named list of global parameters including all arguments as elements
#' @keywords initialize, set up, global, parameter
#' @export

set_g_p <- function(tmax=100, dem.st=FALSE, sdd.st=TRUE, bank=TRUE, n.cores=4, 
                    lc.r=100, lc.c=100, n.lc=6, N.p.t0=10,
                    p.f=c(0.9, 0.1, 0.29, 0.23, 0.2, 0.3),
                    mu=c(200, 100, 40, 20, 20, 10),
                    gamma=2.3, 
                    m=4, 
                    p.c=c(0.3, 0.1, 0.2, 0.2, 0.2, 0.1),
                    sdd.rate=0.1, 
                    sdd.max=15, 
                    bird.hab=c(.35, .35, 0.05, 0.1, 0.1, 0.05), 
                    n.ldd=1,
                    s.c=0.6,
                    s.B=0.3, 
                    s.M=c(0.9, 0.1, 0.6, 0.6, 0.6, 0.6),
                    s.N=c(1, 1, 1, 1, 1, 1),
                    K=c(750, 10, 100, 100, 300, 100),
                    g.D=0.5, 
                    g.B=0.5,
                    p=c(0.07, 0.01, 0.08, 0.02, 0.02, 0.03),
                    edges="wall", method="wt.mn") {
  
  g.p <- list(tmax=tmax, dem.st=dem.st, sdd.st=sdd.st, bank=bank,
              n.cores=n.cores, lc.r=lc.r, lc.c=lc.c, n.lc=n.lc, N.p.t0=N.p.t0, 
              p.f=p.f, mu=mu, gamma=gamma, m=m, 
              p.c=p.c, sdd.rate=sdd.rate, sdd.max=sdd.max, 
              bird.hab=bird.hab, n.ldd=n.ldd, 
              s.c=s.c, s.B=s.B, s.M=s.M, s.N=s.N, K=K, 
              g.D=g.D, g.B=g.B, p=p, 
              edges=edges, method=method)
  
  return(g.p)
}




#' Set buckthorn control treatment parameters
#'
#' Set the buckthorn control treatment parameters for the simulation, allowing
#' for individual elements to be reassigned.
#' @param null_ctrl \code{TRUE} Set control parameters to \code{NULL}?
#' @param t.trt \code{30} Year to start treatments
#' @param add.owners \code{FALSE} Do owners treat every year once starting a
#'   particular treatment?
#' @param grd.i \code{NULL} Vector of cell IDs to receive ground treatments. If
#'   \code{NULL}, then \code{pTrt.grd * ncell} cells are assigned randomly with
#'   the \link{trt_assign} function
#' @param man.i \code{NULL} Vector of cell IDs to receive manual treatments. If
#'   \code{NULL}, then \code{pTrt.man * ncell} cells are assigned randomly with
#'   the \link{trt_assign} function
#' @param chg.i \code{NULL} Vector of cell IDs to receive land cover changes. If
#'   \code{NULL}, then \code{pChg * ncell} cells are assigned randomly with the
#'   \link{trt_assign} function
#' @param pTrt.grd \code{0.05} Proportion of cells with ground treatments in
#'   each time step
#' @param pTrt.man \code{0.05} Proportion of cells with manual treatments in
#'   each time step
#' @param grd.trt \code{Lit=0.005, Cov=0.01, Com=0.00001} Named vector with
#'   ground treatments and associated seedling establishment probabilities
#' @param man.trt \code{c(M=0.1, C=0.3, MC=0.8)} Named vector with manual
#'   treatments and associated mortality (=success) rates
#' @param lc.chg \code{TRUE} Does land cover change across years?
#' @param pChg \code{0.0001} Proportion of cells with land cover change each
#'   year
#' @return Named list of control parameters including all arguments as elements
#'   unless \code{null_ctrl==TRUE}, in which case the function returns
#'   \code{NULL}
#' @keywords initialize, set up, control, treatment, parameter
#' @export

set_control_p <- function(null_ctrl=TRUE, t.trt=30, add.owners=FALSE,
                          grd.i=NULL, man.i=NULL, chg.i=NULL,
                          pTrt.grd=0.05, pTrt.man=0.05,
                          grd.trt=c(Lit=0.005, Cov=0.01, Com=0.00001),
                          man.trt=c(M=0.1, C=0.3, MC=0.8),
                          lc.chg=TRUE, pChg=0.0001) {
  
  if(null_ctrl) {
    control.p <- NULL
  } else {
    control.p <- list(t.trt=t.trt, add.owners=add.owners,
                      grd.i=grd.i, man.i=man.i, chg.i=chg.i,
                      pTrt.grd=pTrt.grd, pTrt.man=pTrt.man,
                      grd.trt=grd.trt, man.trt=man.trt,
                      lc.chg=lc.chg, pChg=pChg)
  }
  
  return(control.p)
}
