##' Iterated Adapted Bagged Filter (IABF)
##'
##' An algorithm for estimating the parameters of a spatiotemporal partially-observed Markov process.
##' Running \code{iabf} causes the algorithm to perform a specified number of iterations of adapted simulations with parameter perturbation and parameter resamplings.
##' At each iteration, adapted simulations are performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' After cycling through the data, each replicate's weight is calculated and is used to rank the bootstrap replictates. The highest ranking replicates are recycled into the next iteration.
##' This extra variability introduced through parameter perturbation effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##'
##' @name iabf
##' @rdname iabf
##' @include spatPomp_class.R abf.R
##' @family particle filter methods
##' @family \pkg{spatPomp} parameter estimation methods
##'
##' @inheritParams pomp::mif2
##' @param Nabf The number of iterations to perform
##' @return
##' Upon successful completion, \code{iabf} returns an object of class
##' \sQuote{iabfd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ extracts the point estimate }
##' }
##'
NULL

rw.sd <- pomp:::safecall


## define the iabfd_spatPomp class
setClass(
  'iabfd_spatPomp',
  contains='abfd_spatPomp',
  slots=c(
    Nabf = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix',
    resamp_frac = 'numeric',
    paramMatrix = 'array'
  )
)

setClass(
  "adapted_replicate_extended",
  contains='adapted_replicate',
  slots=c(
    param = 'numeric',
    log_wp_last = 'numeric',
    log_wm_last = 'numeric'
  ),
  prototype=prototype(
    param=numeric(),
    log_wp_last = numeric(),
    log_wm_last = numeric()
  )
)

h_abf_internal <- function (object, params, Np,
                          abfiter, nbhd, cooling.fn, rw.sd,
                          tol = (1e-18)^17, max.fail = Inf,
                          .indices = integer(0), verbose,
                          .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("h_abf_internal"),": ")
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  abfiter <- as.integer(abfiter)
  Np <- as.integer(Np)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  params <- array(data = params,
                  dim = c(length(params),Np[1L]),
                  dimnames = list(param = names(params), rep = NULL))

  # create array to store weights across time
  log_cond_densities <- array(data = numeric(0), dim=c(nunits,Np[1L],ntimes))
  dimnames(log_cond_densities) <- list(unit = 1:nunits, rep = 1:Np[1L], time = 1:ntimes)
  pompLoad(object,verbose=FALSE)
  for (nt in seq_len(ntimes)) {

    ## perturb parameters
    pmag <- cooling.fn(nt,abfiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')

    tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)

    if (nt == 1L) {
      ## get initial states
      x <- rinit(object,params=tparams)
    }

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        x0=x,
        t0=times[nt],
        times=times[nt+1],
        params=tparams,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    ## determine the weights. returns weights which is a nunits by Np by ntimes array
    log_weights <- tryCatch(
      vec_dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=tparams,
        log=TRUE,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    log_cond_densities[,,nt] <- log_weights[,,1]
    log_resamp_weights <- apply(log_weights[,,1,drop=FALSE], 2, function(x) sum(x))
    max_log_resamp_weights <- max(log_resamp_weights)
    # if any particle's resampling weight is zero replace by tolerance
    if(all(is.infinite(log_resamp_weights))) log_resamp_weights <- rep(log(tol), Np[1L])
    else log_resamp_weights <- log_resamp_weights - max_log_resamp_weights
    resamp_weights <- exp(log_resamp_weights)
    gnsi <- FALSE
    #log.rep.weight <- log.rep.weight + log((1/Np[1L])*sum(resamp_weights))

    xx <- tryCatch(
      .Call(
        "iabf_computations",
        x=X,
        params=params,
        Np=Np[nt+1],
        rw_sd=numeric(0),
        predmean=FALSE,
        predvar=FALSE,
        filtmean=FALSE,
        trackancestry=FALSE,
        onepar=FALSE,
        weights=resamp_weights
      ),
      error = function (e) {
        stop(ep,"iabf resampling computation error: ",conditionMessage(e),call.=FALSE)
      }
    )

    x <- xx$states
    if(nt==ntimes) pred_params <- params
    params <- xx$params

    if (nt == ntimes) {
      if (any(resamp_weights>0)) {
        coef(object, transform=TRUE) <- apply(pred_params,1L,weighted.mean,w=resamp_weights)
      } else {
        warning(ep,"filtering failure at last filter iteration, using unweighted mean for ",
                sQuote("coef"),call.=FALSE)
        coef(object, transform=TRUE) <- apply(pred_params,1L,mean)
      }
    }
  }
  log_loc_comb_pred_weights = array(data = numeric(0), dim=c(nunits,Np[1L], ntimes))
  log_wm_times_wp_avg = array(data = numeric(0), dim = c(nunits, ntimes))
  log_wp_avg = array(data = numeric(0), dim = c(nunits, ntimes))
  for (nt in seq_len(ntimes)){
    for (unit in seq_len(nunits)){
      full_nbhd <- nbhd(object, time = nt, unit = unit)
      log_prod_cond_dens_nt  <- rep(0, Np[1])
      if(length(full_nbhd) > 0) log_prod_cond_dens_not_nt <- matrix(0, Np[1], max(1,nt-min(sapply(full_nbhd,'[[',2))))
      else log_prod_cond_dens_not_nt <- matrix(0,Np[1],0)
      for (neighbor in full_nbhd){
        neighbor_u <- neighbor[1]
        neighbor_n <- neighbor[2]
        if (neighbor_n == nt)
          log_prod_cond_dens_nt  <- log_prod_cond_dens_nt + log_cond_densities[neighbor_u, ,neighbor_n]
        else
          log_prod_cond_dens_not_nt[, nt-neighbor_n] <- log_prod_cond_dens_not_nt[, nt-neighbor_n] + log_cond_densities[neighbor_u, ,neighbor_n]
      }
      log_loc_comb_pred_weights[unit, ,nt]  <- sum(apply(log_prod_cond_dens_not_nt, 2, logmeanexp)) + log_prod_cond_dens_nt
    }
  }
  log_wm_last = log_cond_densities[nunits, ,ntimes]
  log_wp_last = log_loc_comb_pred_weights[nunits, , ntimes]
  params_last = params[,1]
  log_wm_times_wp_avg = apply(log_loc_comb_pred_weights + log_cond_densities, c(1,3), FUN = logmeanexp)
  log_wp_avg = apply(log_loc_comb_pred_weights, c(1,3), FUN = logmeanexp)
  pompUnload(object,verbose=FALSE)
  new(
    "adapted_replicate_extended",
    log_wm_times_wp_avg = log_wm_times_wp_avg,
    log_wm_last = log_wm_last,
    log_wp_avg = log_wp_avg,
    log_wp_last = log_wp_last,
    param=params_last,
    Np=as.integer(Np),
    tol=tol
  )
}

iabf_internal <- function (object, Nrep, nbhd, prop, Nabf, Np, rw.sd,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17, max.fail = Inf, rep_bound = 10000,
                           verbose = FALSE, .ndone = 0L,
                           .indices = vector(mode="list", length = Nrep),
                           .paramMatrix = NULL,
                           .gnsi = TRUE, ...) {

  ep <- paste0("in ",sQuote("iabf"),": ")
  verbose <- as.logical(verbose)
  p_object <- pomp(object,...,verbose=verbose)
  object <- new("spatPomp",p_object,
                unit_covarnames = object@unit_covarnames,
                shared_covarnames = object@shared_covarnames,
                runit_measure = object@runit_measure,
                dunit_measure = object@dunit_measure,
                eunit_measure = object@eunit_measure,
                munit_measure = object@munit_measure,
                vunit_measure = object@vunit_measure,
                unit_names=object@unit_names,
                unitname=object@unitname,
                unit_statenames=object@unit_statenames,
                unit_obsnames = object@unit_obsnames)

  gnsi <- as.logical(.gnsi)
  Np <- c(Np,Np[1L])
  Nabf <- as.integer(Nabf)

  cooling.fn <- pomp:::mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )
  # initialize the parameter for each rep
  .indices <- integer(0)
  if (is.null(.paramMatrix)) {
    rep_param_init <- coef(object)
    start <- coef(object)
  } else {
    rep_param_init <- .paramMatrix
    start <- apply(.paramMatrix,1L,mean)
  }
  rep_param_init <- array(data = rep_param_init,
                          dim = c(length(rep_param_init), Nrep),
                          dimnames = list(param = names(rep_param_init), rep = NULL))

  rw.sd <- pomp:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  traces <- array(dim=c(Nabf+1,length(start)+1),
                    dimnames=list(iteration=seq.int(.ndone,.ndone+Nabf),
                                  variable=c('loglik',names(start))))
  traces[1L,] <- c(NA,start)

  rep_param_init <- partrans(object,rep_param_init,dir="toEst",.gnsi=gnsi)

  # iterate the filtering
  mcopts <- list(set.seed=TRUE)
  ntimes = length(time(object))
  nunits = length(unit_names(object))
  for (n in seq_len(Nabf)) {
    # begin multi-threaded
    mult_rep_output <- list()
    pmag_init <- cooling.fn(1,n)$alpha*rw.sd[,1]*2
    rep_param_init <- .Call('randwalk_perturbation',rep_param_init,pmag_init,PACKAGE = 'pomp')
    mult_rep_output <- foreach::foreach(i=1:Nrep,
                                        .packages=c("pomp","spatPomp"),
                                        .options.multicore=mcopts) %dopar%  {
      spatPomp:::h_abf_internal(
        object=object,
        params=rep_param_init[,i],
        Np=Np,
        nbhd=nbhd,
        abfiter=.ndone+n,
        cooling.fn=cooling.fn,
        rw.sd=rw.sd,
        tol=tol,
        max.fail=max.fail,
        .indices=.indices,
        verbose=verbose
      )
    }
    # end multi-threaded
    # begin single-threaded
    # mult_rep_output <- list()
    # for(i in 1:Nrep){
    #   mult_rep_output[[i]] <- h_abf_internal(
    #     object=object,
    #     params=rep_param_init[,i],
    #     Np=Np,
    #     nbhd=nbhd,
    #     abfiter=.ndone+n,
    #     cooling.fn=cooling.fn,
    #     rw.sd=rw.sd,
    #     tol=tol,
    #     max.fail=max.fail,
    #     .indices=.indices,
    #     verbose=verbose
    #   )
    # }
    # end single threaded
    gnsi <- FALSE
    ## for log-likelihood computation
    cond_loglik <- foreach::foreach(u=seq_len(nunits),
      .combine = 'rbind',
      .packages=c("pomp", "spatPomp"),
      .options.multicore=mcopts) %dopar%
    {
      cond_loglik_u <- array(data = numeric(0), dim=c(ntimes))
      for (n in seq_len(ntimes)){
        log_mp_sum = logmeanexp(vapply(mult_rep_output,
                                       FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,n]),
                                       FUN.VALUE = 1.0))
        log_p_sum = logmeanexp(vapply(mult_rep_output,
                                      FUN = function(rep_output) return(rep_output@log_wp_avg[u,n]),
                                      FUN.VALUE = 1.0))
        cond_loglik_u[n] = log_mp_sum - log_p_sum
      }
      cond_loglik_u
    }
    # for parameter estimation
    rep_loglik_un <- foreach::foreach(u=seq_len(nunits),
                                    .combine = function(...) abind::abind(..., along=3),
                                    .packages=c("pomp", "spatPomp"),
                                    .options.multicore=mcopts) %dopar%
      {
        cond_loglik_u <- array(data = numeric(0), dim=c(ntimes,Nrep))
        for (n in seq_len(ntimes)){
          rep_filt_weight_un_rep = vapply(mult_rep_output,
                                          FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,n]),
                                          FUN.VALUE = 1.0) -
                                   vapply(mult_rep_output,
                                          FUN = function(rep_output) return(rep_output@log_wp_avg[u,n]),
                                          FUN.VALUE = 1.0)
          cond_loglik_u[n,] = rep_filt_weight_un_rep
        }
        cond_loglik_u
      }
    loglik_rep <- apply(rep_loglik_un, MARGIN = 2, FUN = sum)

    ## parameter selection
    max_loglik_rep <- max(loglik_rep)
    # if any particle's resampling weight is zero replace by tolerance
    if(all(is.infinite(loglik_rep))) loglik_rep <- rep(log(tol), Nrep)
    else loglik_rep <- loglik_rep - max_loglik_rep

    param_swarm <- foreach::foreach(i=seq_len(Nrep),
                                    .combine = 'cbind',
                                    .packages=c("pomp", "spatPomp"),
                                    .options.multicore=mcopts) %dopar%
      {
        mult_rep_output[[i]]@param
      }
    top_indices <- which(loglik_rep > quantile(loglik_rep, p = 1-prop[n]))
    new_indices <- rep_len(top_indices, length.out = Nrep)
    rep_param_init <- param_swarm[,new_indices]

    coef(object) <- partrans(object,
                             apply(rep_param_init, 1, mean),
                             dir="fromEst", .gnsi = .gnsi)

    traces[n+1,-c(1)] <- coef(object)
    traces[n+1,c(1)] <- sum(cond_loglik)

    if (n != Nabf){
      rm(param_swarm)
      rm(mult_rep_output)
      rm(top_indices)
      rm(cond_loglik)
      rm(rep_loglik_un)
    }

    if (verbose) cat("iabf iteration",n,"of",Nabf,"completed\n")
  }

  # parameter swarm to be outputted
  param_swarm <- partrans(object,rep_param_init,dir="fromEst",
                                .gnsi=gnsi)

  pompUnload(object,verbose=FALSE)

  new(
    "iabfd_spatPomp",
    object,
    Nabf=as.integer(Nabf),
    Nrep=as.integer(Nrep),
    Np=as.integer(Np),
    resamp_frac=prop,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces,
    paramMatrix=param_swarm,
    loglik=sum(cond_loglik)
  )
}

setGeneric(
  "iabf",
  function (object, ...)
    standardGeneric("iabf")
)

##' @name iabf-iabfd_spatPomp
##' @aliases iabf,iabfd_spatPomp-method
##' @rdname iabf
##' @export
setMethod(
  "iabf",
  signature=signature(object="spatPomp"),
  definition = function (object, Nabf = 1, Nrep, nbhd, prop, Np,
                         rw.sd, cooling.type = c("geometric","hyperbolic"),
                         cooling.fraction.50, tol = (1e-18)^17,
                         max.fail = Inf, rep_bound = 10000,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("iabf"),": ")
    if(missing(Nabf))
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    if (Nabf <= 0)
        stop(ep,sQuote("Nabf")," must be a positive integer",call.=FALSE)
    ntimes <- length(time(object))
    if (missing(Np))
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    else if (!is.numeric(Np))
      stop(ep,sQuote("Np"),
           " must be a number, a vector of numbers, or a function",
           call.=FALSE)
    if (length(Np)==1) {
      Np <- rep(Np,times=ntimes)
    } else if (length(Np)>ntimes) {
      if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
        warning(ep,"Np[k] ignored for k > ntimes",call.=FALSE)
      }
      Np <- head(Np,ntimes)
    }
    if (any(Np <= 0))
      stop(ep,"number of particles, ",
           sQuote("Np"),", must always be positive",call.=FALSE)
    if(missing(Nrep))
      stop(ep,"number of replicates, ",
           sQuote("Nrep"),", must be specified!",call.=FALSE)
    if(missing(prop))
      stop(ep,"top proportion of Nrep (by weight) to be resampled, ",
           sQuote("prop"),", must be specified!",call.=FALSE)
    if(length(prop) == 1) prop <- rep(prop,Nabf)
    if(length(prop) != 1 && length(prop)!= Nabf)
      stop(ep,"top proportion of Nrep (by weight) to be resampled, ",
           sQuote("prop"),", must be of length 1 or ", sQuote("Nabf"), "!",call.=FALSE)


    if (missing(rw.sd))
      stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
    cooling.type <- match.arg(cooling.type)

    cooling.fraction.50 <- as.numeric(cooling.fraction.50)
    if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
      stop(ep,sQuote("cooling.fraction.50"),
           " must be in (0,1]",call.=FALSE)

    iabf_internal(
      object=object,
      Nrep=Nrep,
      prop=prop,
      Nabf=Nabf,
      nbhd=nbhd,
      Np=Np,
      rw.sd=rw.sd,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      max.fail=max.fail,
      rep_bound = rep_bound,
      verbose=verbose,
      ...
    )
  }
)
