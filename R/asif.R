##' Adapted Simulation Island Filter (ASIF)
##'
##' An algorithm for estimating the filter distribution of a spatiotemporal partially-observed Markov process (SpatPOMP for short).
##' Running \code{asif} causes the algorithm to run independent island jobs which each yield an imperfect adapted simulation. Simulating from the "adapted filter"
##' distribution runs into a curse of dimensionality (COD) problem, which is mitigated by keeping particles in each island close to each other through resampling down
##' to one particle per island at each observation time point.
##' The adapted simulations are then weighted in a way that tries to avert COD by making a weak coupling assumption to get an approximate filter distribution.
##' As a by-product, we also get a biased estimate of the likelihood of the data.
##'
##' @name asif
##' @rdname asif
##' @include spatPomp_class.R generics.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @inheritParams spatPomp
##' @inheritParams pomp::pfilter
##' @param object A \code{spatPomp} object.
##' @param params A parameter set for the spatiotemporal POMP. If missing, \code{asif} will attempt to run using \code{coef(object)}
##' @param Np The number of particles used within each island for the adapted simulations.
##' @param nbhd A neighborhood function with three arguments: \code{object}, \code{time} and \code{unit}. The function should return a \code{list} of two-element vectors. The list output of
##' \code{nbhd(u,n)} consists of vectors \code{c(a,b)} where \eqn{(a,b)} is a neighbor of \code{(u,n)} in space-time.
##' @param islands The number of islands for the adapted simulations.
##' @param tol If the resampling weight for a particle is zero due to floating-point precision issues, it is set to the value of \code{tol} since resampling has to be done.
##' @examples
##' # Create a simulation of a BM using default parameter set
##' b <- bm(U=3, N=10)
##'
##' # Create a neighborhood function mapping a point in space-time to a list of ``neighboring points" in space-time
##' bm_nbhd <- function(object, time, unit) {
##'   nbhd_list = list()
##'   if(time > 1 && unit > 1) nbhd_list = c(nbhd_list, list(c(unit - 1, time - 1)))
##'   return(nbhd_list)
##' }
##'
##' # Run ASIF specified number of Monte Carlo islands and particles per island
##' asifd.b <- asif(b, islands = 50, Np = 10, nbhd = bm_nbhd)
##'
##' # Get the likelihood estimate from ASIF
##' logLik(asifd.b)
##'
##' # Compare with the likelihood estimate from Particle Filter
##' pfd.b <- pfilter(b, Np = 500)
##' logLik(pfd.b)
##' @return
##' Upon successful completion, \code{asif} returns an object of class
##' \sQuote{asifd.pomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields a biased estimate of the log-likelihood of the data under the model. }
##' }
##'
NULL

setClass(
  "island.spatPomp",
  slots=c(
    wm.times.wp.avg="array",
    wp.avg="array",
    Np="integer",
    tol="numeric"
  ),
  prototype=prototype(
    wm.times.wp.avg=array(data=numeric(0),dim=c(0,0)),
    wp.avg=array(data=numeric(0),dim=c(0,0)),
    Np=as.integer(NA),
    tol=as.double(NA)
  )
)
setClass(
  "asifd.spatPomp",
  contains="spatPomp",
  slots=c(
    Np="integer",
    tol="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    Np=as.integer(NA),
    tol=as.double(NA),
    loglik=as.double(NA)
  )
)
asif.internal <- function (object, params, Np, nbhd, tol, .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("asif"),": ")
  verbose = FALSE
  if(missing(nbhd))
    stop(ep,sQuote("nbhd")," must be specified for the spatPomp object",call.=FALSE)
  object <- as(object,"spatPomp")
  pompLoad(object,verbose)
  gnsi <- as.logical(.gnsi)

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  if (missing(tol))
    stop(ep,sQuote("tol")," must be specified",call.=FALSE)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(object@units)

  if (missing(Np)) {
    if (is.matrix(params)) {
      Np <- ncol(params)
    } else {
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    }
  }
  if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np")," is a function, ",
             "it must return a single positive integer",call.=FALSE)
      }
    )
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(ep,sQuote("Np")," must have length 1 or length ",ntimes+1,call.=FALSE)
  if (any(Np<=0))
    stop(ep,"number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
  if (!is.numeric(Np))
    stop(ep,sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
  Np <- as.integer(Np)
  if (is.matrix(params)) {
    if (!all(Np==ncol(params)))
      stop(ep,"when ",sQuote("params")," is provided as a matrix, do not specify ",
           sQuote("Np"),"!",call.=FALSE)
  }

  if (NCOL(params)==1) {        # there is only one parameter vector
    one.par <- TRUE
    coef(object) <- params     # set params slot to the parameters
    params <- as.matrix(params)
  }

  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  ## returns an nvars by nsim matrix
  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  loglik <- rep(NA,ntimes)

  # create array to store weights across time
  cond.densities <- array(data = numeric(0), dim=c(nunits,Np[1L],ntimes))
  dimnames(cond.densities) <- list(unit = 1:nunits, rep = 1:Np[1L], time = 1:ntimes)

  for (nt in seq_len(ntimes)) { ## main loop
    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        x0=x,
        t0=times[nt],
        times=times[nt+1],
        params=params,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    ## determine the weights
    weights <- tryCatch(
      vec_dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=params,
        log=FALSE,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    weights[weights == 0] <- tol
    cond.densities[,,nt] <- weights[,,1]
    resamp_weights <- apply(weights[,,1,drop=FALSE], 2, function(x) prod(x))
    # if any particle's resampling weight is zero replace by tolerance
    if(all(resamp_weights == 0)) resamp_weights <- rep(tol, Np[1L])
    gnsi <- FALSE

    ## do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        "asif_computations",
        x=X,
        params=params,
        Np=Np[nt+1],
        trackancestry=FALSE,
        weights=resamp_weights
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )

    x <- xx$states
    params <- xx$params

    if (verbose && (nt%%5==0))
      cat("asif timestep",nt,"of",ntimes,"finished\n")
  } ## end of main loop

  # compute locally combined pred. weights for each time, unit and particle
  loc.comb.pred.weights = array(data = numeric(0), dim=c(nunits,Np[1L], ntimes))
  wm.times.wp.avg = array(data = numeric(0), dim = c(nunits, ntimes))
  wp.avg = array(data = numeric(0), dim = c(nunits, ntimes))
  for (nt in seq_len(ntimes)){
      for (unit in seq_len(nunits)){
          full_nbhd <- nbhd(object, nt, unit)
          prod_cond_dens_nt  <- rep(1, Np[1])
          prod_cond_dens_not_nt <- matrix(1, Np[1], nt-1)
          for (neighbor in full_nbhd){
              neighbor_u <- neighbor[1]
              neighbor_n <- neighbor[2]
              if (neighbor_n == nt)
                  prod_cond_dens_nt  <- prod_cond_dens_nt * cond.densities[neighbor_u, ,neighbor_n]
              else
                  prod_cond_dens_not_nt[, neighbor_n] <- prod_cond_dens_not_nt[, neighbor_n] * cond.densities[neighbor_u, ,neighbor_n]
          }
          loc.comb.pred.weights[unit, ,nt]  <- prod(apply(prod_cond_dens_not_nt, 2, mean))*prod_cond_dens_nt
      }
  }
  wm.times.wp.avg = apply(loc.comb.pred.weights * cond.densities, c(1,3), FUN = mean)
  wp.avg = apply(loc.comb.pred.weights, c(1,3), FUN = mean)

  pompUnload(object,verbose=verbose)
  new(
    "island.spatPomp",
    wm.times.wp.avg = wm.times.wp.avg,
    wp.avg = wp.avg,
    Np=as.integer(Np),
    tol=tol
  )


}
##' @name asif-spatPomp
##' @aliases asif,spatPomp-method
##' @rdname asif
##' @export
setMethod(
  "asif",
  signature=signature(object="spatPomp"),
  function (object, islands, Np, nbhd, params,
           tol = (1e-18)^9,
           ...) {
   if (missing(params)) params <- coef(object)
   ## single thread for testing
   # single_island_output <- asif.internal(
   #  object=object,
   #  params=params,
   #  Np=Np,
   #  nbhd = nbhd,
   #  tol=tol,
   #  ...)
   # return(single_island_output)
   ## end single thread for testing
   ## cores <- parallel:::detectCores() - 1
   #
   ## foreach now registered outside asif
   ## doParallel::registerDoParallel(cores = NULL)
   #
   ## begin multi-thread code
   mcopts <- list(set.seed=TRUE)
   # set.seed(396658101,kind="L'Ecuyer")
   mult_island_output <- foreach::foreach(i=1:islands,
       .packages=c("pomp","spatPomp"),
       .options.multicore=mcopts) %dopar%  spatPomp:::asif.internal(
     object=object,
     params=params,
     Np=Np,
     nbhd = nbhd,
     tol=tol,
     ...
     )
   ntimes = length(time(object))
   nunits = length(spat_units(object))
   # compute sum (over all islands) of w_{d,n,i}^{P} for each (d,n)
   island_mp_sums = array(data = numeric(0), dim = c(nunits,ntimes))
   island_p_sums = array(data = numeric(0), dim = c(nunits, ntimes))
   cond.loglik = array(data = numeric(0), dim=c(nunits, ntimes))
   for (i in seq_len(nunits)){
    for (j in seq_len(ntimes)){
      mp_sum = 0
      p_sum = 0
      for (k in seq_len(islands)){
        mp_sum = mp_sum + mult_island_output[[k]]@wm.times.wp.avg[i,j]
        p_sum = p_sum + mult_island_output[[k]]@wp.avg[i,j]
      }
      cond.loglik[i,j] = log(mp_sum) - log(p_sum)
    }
   }
   # end multi-threaded code
   new(
      "asifd.spatPomp",
      object,
      Np=as.integer(Np),
      tol=tol,
      loglik=sum(cond.loglik)
     )
  }
)

