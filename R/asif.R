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
##' @include spatpomp_class.R generics.R
##' @family particle filter methods
##' @family \pkg{spatpomp} filtering methods
##'
##'
##' @inheritParams spatpomp
##' @inheritParams pfilter
##' @param object A \code{spatpomp} object.
##' @param params A parameter set for the spatiotemporal POMP.
##' @param Np The number of particles used within each island for the adapted simulations.
##' @param nbhd A function which when given a point in space-time \code{(d,n)} will output a 0-1 matrix of points in space-time considered to be
##' in the neighborhood of \code{(d,n)}.
##' @param islands The number of islands for the adapted simulations.
##' @return
##' Upon successful completion, \code{asif} returns an object of class
##' \sQuote{asifd.pomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{loglik}}}{ yields a biased estimate of the log-likelihood of the data under the model. }
##' }
##'
NULL

setClass(
  "island.spatpomp",
  contains="spatpomp",
  slots=c(
    loc.comb.pred.weights="array",
    cond.densities="array",
    #pred.mean="array",
    #pred.var="array",
    #filter.mean="array",
    #filter.traj="array",
    #paramMatrix="array",
    #indices="vector",
    #eff.sample.size="numeric",
    #cond.loglik="numeric",
    #saved.states="list",
    #saved.params="list",
    Np="integer",
    tol="numeric"
    #nfail="integer",
    #loglik="numeric"
  ),
  prototype=prototype(
    loc.comb.pred.weights=array(data=numeric(0),dim=c(0,0)),
    cond.densities=array(data=numeric(0),dim=c(0,0,0)),
    # pred.mean=array(data=numeric(0),dim=c(0,0)),
    # pred.var=array(data=numeric(0),dim=c(0,0)),
    # filter.mean=array(data=numeric(0),dim=c(0,0)),
    # filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    # paramMatrix=array(data=numeric(0),dim=c(0,0)),
    # indices=integer(0),
    # eff.sample.size=numeric(0),
    # cond.loglik=numeric(0),
    # saved.states=list(),
    # saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA)
    #nfail=as.integer(NA),
    #loglik=as.double(NA)
  )
)
setClass(
  "asifd.spatpomp",
  contains="spatpomp",
  slots=c(
    #loc.comb.filter.weights="array",
    #pred.mean="list",
    #pred.var="list",
    #filter.mean="list",
    #filter.traj="list",
    #paramMatrix="array",
    #indices="vector",
    #eff.sample.size="list",
    #cond.loglik="array",
    #saved.states="list",
    #saved.params="list",
    Np="integer",
    tol="numeric",
    #nfail="list",
    loglik="numeric"
  ),
  prototype=prototype(
    #loc.comb.filter.weights=array(data=numeric(0),dim=c(0,0,0,0)),
    #pred.mean=list(),
    #pred.var=list(),
    #filter.mean=list(),
    #filter.traj=list(),
    #paramMatrix=array(data=numeric(0),dim=c(0,0)),
    #indices=integer(0),
    #eff.sample.size=list(),
    #cond.loglik=array(data=numeric(0),dim=c(0,0)),
    #saved.states=list(),
    #saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    #nfail=list(),
    loglik=as.double(NA)
  )
)
asif.internal <- function (object, params, Np,
                              nbhd,
                              tol, max.fail,
                              pred.mean = FALSE,
                              pred.var = FALSE,
                              filter.mean = FALSE,
                              filter.traj = FALSE,
                              cooling, cooling.m,
                              verbose = FALSE,
                              save.states = FALSE,
                              save.params = FALSE,
                              .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("asif"),": ")
  if(missing(nbhd))
    stop(ep,sQuote("nbhd")," must be specified for the spatpomp object",call.=FALSE)
  ep <- paste0("in ",sQuote("asif"),": ")
  object <- as(object,"spatpomp")
  pompLoad(object,verbose=verbose)
  gnsi <- as.logical(.gnsi)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  verbose <- as.logical(verbose)
  save.states <- as.logical(save.states)
  save.params <- as.logical(save.params)

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  if (missing(tol))
    stop(ep,sQuote("tol")," must be specified",call.=FALSE)

  one.par <- FALSE
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
  #init.x <- init.state(object,params=params,nsim=Np[1L])
  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states | filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (save.params) {
    pparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  } else {
    pparticles <- list()
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(
      data=0,
      dim=c(nvars,1,ntimes+1),
      dimnames=list(
        variable=statenames,
        rep=1,
        time=times)
    )
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

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

    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
      if (length(problem.indices)>0) {  # state variables
        stop(
          ep,"non-finite state variable(s): ",
          paste(rownames(X)[problem.indices],collapse=', '),
          call.=FALSE
        )
      }
    }
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
    # if any particle's resampling weight is zero divide out it's weight vector by the smallest component
    if(all(resamp_weights == 0)) resamp_weights <- rep(tol, Np[1L])
    gnsi <- FALSE
    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- tryCatch(
      .Call(
        "asif_computations",
        x=X,
        params=params,
        Np=Np[nt+1],
        rw_sd=numeric(0),
        predmean=pred.mean,
        predvar=pred.var,
        filtmean=filter.mean,
        trackancestry=filter.traj,
        onepar=one.par,
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

  # compute locally combined pred. weights for each time and unit
  loc.comb.pred.weights = array(data = numeric(0), dim=c(nunits,ntimes))
  for (nt in seq_len(ntimes)){
    for (unit in seq_len(nunits)){
      prod_over_times = 1
      full_nbhd = nbhd(object, nt, unit)
      for (prev_t in 1:nt){
        if(prev_t == nt){
          if(unit == 1) next
          time_sum = 0
          for(pp in seq_len(Np[1])){
            part_prod = 1
            for (prev_u in 1:unit-1){
              if (prev_u == 0) next
              if (full_nbhd[prev_u, prev_t]){
                part_prod = part_prod*cond.densities[prev_u, pp, prev_t]
              }
            }
            time_sum = time_sum + part_prod
          }
          time_avg = time_sum/Np[1]
          prod_over_times = prod_over_times * time_avg
        }
        else{
          time_sum = 0
          for(pp in seq_len(Np[1])){
            part_prod = 1
            for (prev_u in 1:nunits){
              if (full_nbhd[prev_u, prev_t]){
                part_prod = part_prod*cond.densities[prev_u, pp, prev_t]
              }
            }
            time_sum = time_sum + part_prod
          }
          time_avg = time_sum/Np[1]
          prod_over_times = prod_over_times * time_avg
        }
      }
      loc.comb.pred.weights[unit,nt] = prod_over_times
    }
  }
  pompUnload(object,verbose=verbose)
  new(
    "island.spatpomp",
    object,
    loc.comb.pred.weights = loc.comb.pred.weights,
    cond.densities = cond.densities,
    #pred.mean=pred.m,
    #pred.var=pred.v,
    #filter.mean=filt.m,
    #filter.traj=filt.t,
    #paramMatrix=array(data=numeric(0),dim=c(0,0)),
    #eff.sample.size=eff.sample.size,
    #cond.loglik=loglik,
    #saved.states=xparticles,
    #saved.params=pparticles,
    Np=as.integer(Np),
    tol=tol
    #nfail=as.integer(nfail),
    #loglik=sum(loglik)
  )


}
##' @name asif-spatpomp
##' @aliases asif,spatpomp-method
##' @rdname asif
##' @export
setMethod(
  "asif",
  signature=signature(object="spatpomp"),
  function (object, params, Np, nbhd, islands,
           tol = (1e-18)^9,
           max.fail = Inf,
           pred.mean = FALSE,
           pred.var = FALSE,
           filter.mean = FALSE,
           filter.traj = FALSE,
           save.states = FALSE,
           save.params = FALSE,
           verbose = getOption("verbose"),
           ...) {
   if (missing(params)) params <- coef(object)
   ## single thread for testing
   # single_island_output <- asif.internal(
   #  object=object,
   #  params=params,
   #  Np=Np,
   #  nbhd = nbhd,
   #  tol=tol,
   #  max.fail=max.fail,
   #  pred.mean=pred.mean,
   #  pred.var=pred.var,
   #  filter.mean=filter.mean,
   #  filter.traj=filter.traj,
   #  save.states=save.states,
   #  save.params=save.params,
   #  verbose=verbose,
   #  ...)
   # return(single_island_output)
   ## end single thread for testing
   ## begin multi-thread code
   require(doParallel)
   require(parallel)
   # cores <- parallel:::detectCores() - 1
   registerDoParallel(cores = NULL)
   mcopts <- list(set.seed=TRUE)
   # set.seed(396658101,kind="L'Ecuyer")
   mult_island_output <- foreach(i=1:islands, .options.multicore=mcopts) %dopar%  asif.internal(
     object=object,
     params=params,
     Np=Np,
     nbhd = nbhd,
     tol=tol,
     max.fail=max.fail,
     pred.mean=pred.mean,
     pred.var=pred.var,
     filter.mean=filter.mean,
     filter.traj=filter.traj,
     save.states=save.states,
     save.params=save.params,
     verbose=verbose,
     ...
     )
   ntimes = length(time(object))
   nunits = length(unit(object))
   # compute sum (over all islands) of w_{d,n,i}^{P} for each (d,n)
   island_weight_sums = array(data = numeric(0), dim = c(nunits,ntimes))
   for (i in seq_len(nunits)){
    for (j in seq_len(ntimes)){
      weight_sum = 0
      for (k in seq_len(islands)){
        # weight_sum = weight_sum + mult_island_output[[k]][[2]][i,j]
        weight_sum = weight_sum + mult_island_output[[k]]@loc.comb.pred.weights[i,j]
      }
      island_weight_sums[i,j] = weight_sum
    }
   }
   loc.comb.filter.weights = array(data = numeric(0), dim=c(nunits, ntimes, islands, Np))
   marg.filter.probs = array(data = numeric(0), dim=c(nunits, ntimes, islands, Np))
   # compute w_{d,n,i,j}^LCF
   for(i in seq_len(nunits)){
    for(j in seq_len(ntimes)){
      for(k in seq_len(islands)){
        for(l in seq_len(Np)){
           loc.comb.filter.weights[i,j,k,l] = mult_island_output[[k]]@cond.densities[i,l,j] * (mult_island_output[[k]]@loc.comb.pred.weights[i,j]) / island_weight_sums[i,j]
        }
      }
    }
   }
   # end multi-threaded code
   #
   # compute conditional log-likelihood estimate
   cond.loglik = array(data = numeric(0), dim=c(nunits, ntimes))
   for(i in seq_len(nunits)){
    for(j in seq_len(ntimes)){
      space_time_sum = 0
      for(k in seq_len(islands)){
        for(l in seq_len(Np)){
           space_time_sum = space_time_sum + loc.comb.filter.weights[i,j,k,l]
        }
      }
      cond.loglik[i,j] = log((1/Np)*space_time_sum)
    }
   }

   new(
      "asifd.spatpomp",
      object,
      #loc.comb.filter.weights = loc.comb.filter.weights,
      #pred.mean=lapply(mult_island_output, "slot", "pred.mean"),
      #pred.var=lapply(mult_island_output, "slot", "pred.var"),
      #filter.mean=lapply(mult_island_output, "slot", "filter.mean"),
      #filter.traj=lapply(mult_island_output, "slot", "filter.traj"),
      #paramMatrix=array(data=numeric(0),dim=c(0,0)),
      #eff.sample.size=lapply(mult_island_output, "slot", "eff.sample.size"),
      #cond.loglik=cond.loglik,
      #saved.states=lapply(mult_island_output, "slot", "saved.states"),
      #saved.params=lapply(mult_island_output, "slot", "saved.params"),
      Np=as.integer(Np),
      tol=tol,
      #nfail=lapply(mult_island_output, "slot", "nfail"),
      loglik=sum(cond.loglik)
     )
  }
)

