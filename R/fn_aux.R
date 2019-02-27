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



#' Get the id of the cell containing specified coordinates
#' 
#' @param lc.df Dataframe or tibble with columns 'lat' and 'lon'
#' @param pt_coord Vector with longitude and latitude of a point
#' @return id from lc.df of cell containing the point
#' @export

get_pt_id <- function(lc.df, pt_coord) {
  cell_side <- mean(diff(sort(unique(lc.df$lon))))
  lc.df$id[which(abs(lc.df$lon-pt_coord[1]) < cell_side/2 & 
                   abs(lc.df$lat-pt_coord[2]) < cell_side/2)]
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




#' Expand land cover-based parameter ranges (deprecated)
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
  library(tidyverse)
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




#' Expand to all combinations of canopy parameter ranges (deprecated)
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




#' Generate landscape grid
#'
#' Create rectangular landscape grid with x and y coordinates, inbound
#' indicator, and land cover composition
#' @param in.file File (.csv) containing land cover composition and coordinates
#' @param lon \code{"lon"} Column name containing x (longitude) coordinates
#' @param lat \code{"lat"} Column name containing x (latitude) coordinates
#' @param col.inc Indexes of land cover (or other covariate) columns
#' @param out.file \code{NULL} File to store gridded output
#' @return Dataframe with length(x)*length(y) rows and columns \code{x, y, x_y,
#'   inbd} and land cover or covariates
#' @keywords landscape, grid, initialize
#' @export

make_grid <- function(in.file, x.="lon", y.="lat", col.inc, out.file=NULL) {
  library(tidyverse)
  raw.df <- read_csv(in.file) %>%
    mutate(x=as.numeric(as.factor(.[[x.]])),
           y=as.numeric(factor(.[[y.]], levels=rev(levels(factor(.[[y.]]))))),
           x_y=paste(x, y, sep="_"))
  lc.rct <- as.tibble(expand.grid(x=1:max(raw.df$x),
                                  y=1:max(raw.df$y))) %>%
    mutate(x_y=paste(x, y, sep="_"))
  match.order <- match(lc.rct$x_y, raw.df$x_y)
  for(i in col.inc) {
    lc.rct[[names(raw.df)[i]]] <- raw.df[[i]][match.order]
    lc.rct[[names(raw.df)[i]]][is.na(match.order)] <- 0
  }
  lc.rct$lon <- raw.df[[x.]][match.order]
  lc.rct$lat <- raw.df[[y.]][match.order]
  lc.rct$inbd <- !is.na(match.order)
  if(!is.null(out.file)) saveRDS(lc.rct, paste0(out.file))
  return(lc.rct)
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
#' @param s.N Vector \code{length=n.lc} of annual adult survival rates or vector
#'   of slopes corresponding with columns in lc.df
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
#' @param p.trt_OpnOnly \code{FALSE} Do ground cover treatments only apply to
#'   Open land cover types?
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
                   p.trt=NULL, edges="wall", method="wt.mn", p.trt_OpnOnly=F) {
  
  library(tidyverse)
  
  if(method=="wt.mn") {
    lc.mx <- as.matrix(select(lc.df, 
                        -one_of("x", "y", "x_y", "inbd", "id", "id.in", "lat", "lon")))
    nLC <- ncol(lc.mx)
    # scalar = same value for all LC categories
    if(length(K)==1) K <- rep(K, nLC)
    if(length(s.N)==1) s.N <- rep(s.N, nLC)
    if(length(s.M)==1) s.M <- rep(s.M, nLC)
    if(length(mu)==1) mu <- rep(mu, nLC)
    if(length(p.f)==1) p.f <- rep(p.f, nLC)
    if(length(p.c)==1) p.c <- rep(p.c, nLC)
    if(length(p)==1) p <- rep(p, nLC)
    # take weighted mean
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
    if(p.trt_OpnOnly) {
      if(method=="wt.mn") {
        for(i in seq_along(p.trt$id)) {
          p.E[p.trt$id[i],] <- lc.mx[p.trt$id[i],] %*% c(p.trt$p[i], p[-1])
        }
      } else if(method=="lm") {
        for(i in seq_along(p.trt$id)) {
          p.E[p.trt$id[i],] <- antilogit(lc.mx[p.trt$id[i],1:length(p)] %*% 
                                           c(p.trt$p[i], p[-1]))
        }
      }
    } else {
      p.E[p.trt$id,] <- p.trt$p
    }
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
#' @param p.0 \code{NULL} Cell IDs with populations at t=0
#' @param N.0 \code{NULL} Initial adult abundance in cell
#' @return Matrix or array of initial abundances with \code{dim=c(ngrid, (n.lc),
#'   m.max)}
#' @keywords initialize, set up
#' @export

pop_init <- function(ngrid, g.p, lc.df, p.0=NULL, N.0=NULL) {
  
  if(is.null(p.0)) p.0 <- sample(lc.df$id[lc.df$inbd], g.p$N.p.t0)
  if(is.null(N.0)) N.0 <- g.p$K/2
  m.max <- max(g.p$m)  # adult age bin
  
  if(length(unique(g.p$m)) == 1) {
    N.init <- matrix(0, ngrid, m.max)  # column for each age class
    N.init[p.0,m.max] <- round(as.matrix(lc.df[lc.df$id %in% p.0, 4:9]) %*% 
                                N.0)
    N.init[p.0,-m.max] <- round(N.init[p.0,m.max]/5)
    
  } else {
    N.init <- array(0, dim=c(ngrid, g.p$n.lc, m.max))
    N.init[p.0,,m.max] <- round(t(t(as.matrix(lc.df[lc.df$id %in% p.0, 4:9])) * 
                                   N.0))
    N.init[p.0,,-m.max] <- round(N.init[p.0,,m.max]/5)
  }
  
  return(N.init)
}




#' Set global parameters
#'
#' Set the global parameters for the simulation, allowing for individual
#' elements to be reassigned. Default parameters assume 8.1 hectare cells.
#' @param tmax \code{100} Number of time steps per simulation
#' @param dem.st \code{FALSE} Include stochasticity in demography?
#' @param sdd.st \code{TRUE} Include stochasticity in short distance dispersal?
#' @param bank \code{TRUE} Include seedbank?
#' @param n.cores \code{4} Number of cores for parallelizing sdd.pr calculation
#' @param lc.r \code{Inf} Maximum number of rows (\code{y}) in landscape
#' @param lc.c \code{Inf} Maximum number of columns (\code{x}) in landscape
#' @param n.lc \code{6} Number of land cover categories
#' @param N.p.t0 \code{10} Number of cells with buckthorn at t=1
#' @param p.f \code{c(0.45, 0, 0.291, 0.309, 0.309, 0.272)} Vector
#'   \code{length=n.lc} of fruiting probabilities
#' @param mu \code{c(1948, 0, 14, 41, 41, 21)} Vector \code{length=n.lc} of mean
#'   fruit per adult
#' @param gamma \code{2.48} Scalar: mean number of seeds per fruit
#' @param m \code{c(3, 3, 7, 7, 7, 7)} Vector \code{length=n.lc} or scalar of
#'   age at first fruiting. Individuals at this age are considered adults
#' @param p.c \code{c(0.165, 0.165, 0.296, 0.252, 0.252, 0.296)} Vector
#'   \code{length=n.lc} of proportion of fruits eaten by birds, with
#'   \code{1-p.c} assumed to drop directly below buckthorn individuals
#' @param sdd.rate \code{0.03775} 1/mn for exponential dispersal kernel (units =
#'   cells); value assumes 20ac (8.1 ha) grid cells
#' @param sdd.max \code{24} Maximum dispersal distance in cells (units = cells);
#'   value assumes 20ac (8.1 ha) grid cells
#' @param bird.hab \code{c(0.32, 0.36, 0.05, 0.09, 0.09, 0.09)} Vector
#'   \code{length=n.lc} of relative bird habitat preferences
#' @param n.ldd \code{19} Number of long distance dispersal events per year
#' @param s.c \code{0.585} Seed viability post-digestion
#' @param s.B \code{0.72} Probability of annual survival in seed bank
#' @param s.M \code{c(0.9, 0. 0.6, 0.6, 0.6, 0.6)} Vector \code{length=n.lc} of
#'   annual juvenile survival rates
#' @param s.N \code{c(1, 1, 1, 1, 1, 1)} Vector \code{length=n.lc} of annual
#'   adult survival rates
#' @param K \code{c(28205, 0, 4162, 4162, 4162, 4162)} Vector (length=n.lc) of
#'   carrying capacities for adults; values assume 20ac (8.1 ha) grid cells
#' @param g.D \code{0} Probability of direct germination (i.e., a seed
#'   germinates in the same year it is produced)
#' @param g.B \code{0.2} Probability of germinating from the seed bank
#' @param p \code{c(0.295, 0, 0.421, 0.082, 0.082, 0.23)} Vector
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
                    lc.r=Inf, lc.c=Inf, n.lc=6, N.p.t0=10,
                    p.f=c(0.45, 0, 0.291, 0.309, 0.309, 0.272),
                    mu=c(1948, 0, 14, 41, 41, 21),
                    gamma=2.48, 
                    m=c(3, 3, 7, 7, 7, 7), 
                    p.c=c(0.165, 0.165, 0.296, 0.252, 0.252, 0.296),
                    sdd.rate=0.03775, 
                    sdd.max=24, 
                    bird.hab=c(0.32, 0.36, 0.05, 0.09, 0.09, 0.09), 
                    n.ldd=19,
                    s.c=0.585,
                    s.B=0.72, 
                    s.M=c(0.9, 0, 0.6, 0.6, 0.6, 0.6),
                    s.N=c(1, 1, 1, 1, 1, 1),
                    K=c(28205, 0, 4162, 4162, 4162, 4162),
                    g.D=0, 
                    g.B=0.2,
                    p=c(0.295, 0, 0.421, 0.082, 0.082, 0.23),
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




#' Set buckthorn control (management) treatment parameters
#'
#' Set the buckthorn control treatment parameters for the simulation, allowing
#' for individual elements to be reassigned. Treatments fall into two broad
#' categories: ground cover and manual. Ground cover treatments reduce seedling
#' recruitment (i.e., establishment probabilities, \code{set_g_p()$p}) by
#' covering the ground with litter or cover crops, or by compacting the soil.
#' When these treatments are enacted in a cell, the establishment probability is
#' replaced by the corresponding probability assigned in
#' \code{set_control_p()$grd.trt}. Manual treatments reduce juvenile and adult
#' abundances by mechanical and/or chemical methods. The success rates (i.e.,
#' mortality rates) are assigned in \code{set_control_p()$man.trt}, and when a
#' manual treatment is enacted in a cell, the adult and juvenile abundances in
#' that cell are reduced by the corresponding rate.
#' @param null_ctrl \code{TRUE} Set all control parameters to \code{NULL}?
#' @param t.trt \code{1} Year to start treatments
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
#'   each time step; only used if \code{grd.i=NULL}
#' @param pTrt.man \code{0.05} Proportion of cells with manual treatments in
#'   each time step; only used if \code{man.i=NULL}
#' @param grd.trt \code{c(Lit=0.08, Cov=0.175, Com=0.246)} Named vector with
#'   ground treatments and associated seedling establishment probabilities;
#'   treatments include litter (Lit), cover crops (Cov), and soil compaction
#'   (Com)
#' @param man.trt \code{c(M=0.68, C=0.9, MC=0.97, N=0)} Named vector with manual
#'   treatments and associated mortality (=success) rates; treatments include
#'   mechanical (M), chemical (C), both (MC), or none (N)
#' @param lc.chg \code{FALSE} Does land cover change (i.e., timber harvest)
#'   across years?
#' @param pChg \code{0.0001} Proportion of cells with land cover change each
#'   year
#' @return Named list of control parameters including all arguments as elements
#'   unless \code{null_ctrl==TRUE}, in which case the function returns
#'   \code{NULL}
#' @keywords initialize, set up, control, treatment, parameter
#' @export

set_control_p <- function(null_ctrl=TRUE, t.trt=1, add.owners=FALSE,
                          grd.i=NULL, man.i=NULL, chg.i=NULL,
                          pTrt.grd=0.05, pTrt.man=0.05,
                          grd.trt=c(Lit=0.08, Cov=0.175, Com=0.246),
                          man.trt=c(M=0.68, C=0.9, MC=0.97, N=0),
                          lc.chg=FALSE, pChg=0.0001) {
  
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







