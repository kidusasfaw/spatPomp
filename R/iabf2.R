##' Iterated Adapted Bagged Filter 2 (IABF2)
##'
##' An algorithm for estimating the parameters of a spatiotemporal partially-observed Markov process.
##' Running \code{iabf} causes the algorithm to perform a specified number of iterations of adapted simulations with parameter perturbation and parameter resamplings.
##' At each iteration, adapted simulations are performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' After cycling through the data, each replicate's weight is calculated and is used to rank the bootstrap replictates. The highest ranking replicates are recycled into the next iteration.
##' This extra variability introduced through parameter perturbation effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##'
##' @name iabf2
##' @rdname iabf2
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
  'iabf2d_spatPomp',
  contains='abfd_spatPomp',
  slots=c(
    Nabf = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix',
    paramMatrix = 'array'
  )
)

setClass(
  "adapted_replicate_extended2",
  contains='adapted_replicate',
  slots=c(
    param = 'numeric',
    state = 'numeric',
    prev_weights = 'array'
  ),
  prototype=prototype(
    param=numeric(0),
    state=numeric(0),
    prev_weights = array(data=numeric(0),dim=c(0,0,0))
  )
)

h_abf_internal2 <- function (object,
                             params,
                             states,
                             obs_num,
                             prev_meas_weights,
                             Np,
                             nbhd,
                             abfiter,
                             cooling.fn,
                             rw.sd,
                             tol = (1e-18)^17, max.fail = Inf,
                             .indices = integer(0), verbose,
                             .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("h_abf_internal2"),": ")
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  abfiter <- as.integer(abfiter)
  Np <- as.integer(Np)

  nt <- obs_num
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))

  params <- array(data = params,
                  dim = c(length(params),Np[1L]),
                  dimnames = list(param = names(params), rep = NULL))

  if(!is.null(prev_meas_weights))
    prev_meas_weights <- array(data = prev_meas_weights,
                               dim = c(dim(prev_meas_weights)[1],
                                       dim(prev_meas_weights)[2],
                                       dim(prev_meas_weights)[3]),
                               )

  # create array to store weights across time
  log_cond_densities <- array(data = numeric(0), dim=c(nunits,Np[1L]))
  dimnames(log_cond_densities) <- list(unit = 1:nunits, rep = 1:Np[1L])
  pompLoad(object,verbose=FALSE)

  ## perturb parameters
  pmag <- cooling.fn(nt,abfiter)$alpha*rw.sd[,nt]
  params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')

  tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)

  if (nt == 1L) {
    ## get initial states
    x <- rinit(object,params=tparams)
  } else{
    x <- states
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

  ## determine the weights. returns weights which is a nunits by Np array
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
  log_cond_densities <- log_weights[,,1]
  log_resamp_weights <- apply(log_weights[,,1,drop=FALSE], 2, function(x) sum(x))
  max_log_resamp_weights <- max(log_resamp_weights)
  # if any particle's resampling weight is zero replace by tolerance
  if(all(is.infinite(log_resamp_weights))) log_resamp_weights <- rep(log(tol), Np[1L])
  else log_resamp_weights <- log_resamp_weights - max_log_resamp_weights
  resamp_weights <- exp(log_resamp_weights)
  gnsi <- FALSE

  xx <- tryCatch(
    .Call(
      "iabf_computations",
      x=X,
      params=params,
      Np=Np[1L],
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
  # if(nt==ntimes) pred_params <- params
  params <- xx$params

  # if (nt == ntimes) {
  #   if (any(resamp_weights>0)) {
  #     coef(object, transform=TRUE) <- apply(pred_params,1L,weighted.mean,w=resamp_weights)
  #   } else {
  #     warning(ep,"filtering failure at last filter iteration, using unweighted mean for ",
  #             sQuote("coef"),call.=FALSE)
  #     coef(object, transform=TRUE) <- apply(pred_params,1L,mean)
  #   }
  # }

  log_loc_comb_pred_weights = array(data = numeric(0), dim=c(nunits,Np[1L]))
  log_wm_times_wp_avg = array(data = numeric(0), dim = c(nunits))
  log_wp_avg = array(data = numeric(0), dim = c(nunits))
  for (unit in seq_len(nunits)){
    full_nbhd <- nbhd(object, time = nt, unit = unit)
    log_prod_cond_dens_nt  <- rep(0, Np[1])
    if(length(full_nbhd) > 0) log_prod_cond_dens_not_nt <- matrix(0, Np[1], max(1,nt-min(sapply(full_nbhd,'[[',2))))
    else log_prod_cond_dens_not_nt <- matrix(0,Np[1],0)
    for (neighbor in full_nbhd){
      neighbor_u <- neighbor[1]
      neighbor_n <- neighbor[2]
      if (neighbor_n == nt)
        log_prod_cond_dens_nt  <- log_prod_cond_dens_nt + log_cond_densities[neighbor_u, ]
      else{
        if(!is.null(prev_meas_weights))
          log_prod_cond_dens_not_nt[, nt-neighbor_n] <- log_prod_cond_dens_not_nt[, nt-neighbor_n] +
            prev_meas_weights[neighbor_u,dim(prev_meas_weights)[2]-(nt-neighbor_n-1),]
      }
    }
    log_loc_comb_pred_weights[unit, ]  <- sum(apply(log_prod_cond_dens_not_nt, 2, logmeanexp)) + log_prod_cond_dens_nt
  }
  if(!is.null(prev_meas_weights))
    prev_meas_weights <- abind::abind(prev_meas_weights, log_cond_densities, along = 2)
  else
    prev_meas_weights <- array(log_cond_densities, dim=c(nunits, 1, Np[1]))

  params_last = params[,1]
  log_wm_times_wp_avg = apply(log_loc_comb_pred_weights + log_cond_densities, c(1), FUN = logmeanexp)
  log_wp_avg = apply(log_loc_comb_pred_weights, c(1), FUN = logmeanexp)
  pompUnload(object,verbose=FALSE)
  new(
    "adapted_replicate_extended2",
    log_wm_times_wp_avg = array(log_wm_times_wp_avg, dim=c(nunits,1)),
    log_wp_avg = array(log_wp_avg, dim=c(nunits,1)),
    param=params_last,
    state=x[,1],
    prev_weights=prev_meas_weights,
    Np=as.integer(Np),
    tol=tol
  )
}

iabf_internal2 <- function (object, Nrep, nbhd, Nabf, Np, rw.sd,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17, max.fail = Inf,
                           verbose = FALSE, .ndone = 0L,
                           .indices = integer(0),
                           .paramMatrix = NULL,
                           .gnsi = TRUE, ...) {

  ep <- paste0("in ",sQuote("iabf2"),": ")
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

  rw.sd <- pomp:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  traces <- array(dim=c(Nabf+1,length(start)+1),
                    dimnames=list(iteration=seq.int(.ndone,.ndone+Nabf),
                                  variable=c('loglik',names(start))))
  traces[1L,] <- c(NA,start)

  rep_param_init <- partrans(object,rep_param_init,dir="toEst",.gnsi=gnsi)

  # iterate the filtering
  mcopts <- list(set.seed=TRUE)
  times = time(object)
  ntimes = length(time(object))
  nunits = length(unit_names(object))
  for (n in seq_len(Nabf)) {
    # begin single-threaded
    prev_weights <- NULL
    cond_loglik <- array(dim=c(nunits,ntimes))
    states <- NULL
    for(nt in seq_len(ntimes)){
      mult_rep_output <- list()
      for(i in seq_len(Nrep)){
        mult_rep_output <- c(mult_rep_output, spatPomp:::h_abf_internal2(
          object=object,
          params=rep_param_init,
          states=(if(is.null(states)) NULL else states[,i]),
          obs_num=nt,
          prev_meas_weights=(if(is.null(prev_weights)) NULL else prev_weights[,,,i,drop=FALSE]),
          Np=Np,
          nbhd=nbhd,
          abfiter=.ndone+n,
          cooling.fn=cooling.fn,
          rw.sd=rw.sd,
          tol=tol,
          max.fail=max.fail,
          .indices=.indices,
          verbose=verbose
        ))
      }
      # for the next observation time, how far back do we need
      # to provide the conditional densities, f_{Y_{u,n}|X_{u,n}}?
      max_lookback <- 0
      for(u in seq_len(nunits)){
        farthest_back <- min(sapply(nbhd(unit=u,time=(nt+1)),'[[',2))
        if(nt+1 - farthest_back > max_lookback) max_lookback <- nt+1 - farthest_back
      }
      # update prev_weights for the next observation time
      # THIS NEEDS TIGHTENING UP FOR DIFFERENT KINDS OF NEIGHBORHOODS
      prev_prev_weights <- foreach::foreach(
        i=seq_len(Nrep),
        .combine = function(...) abind::abind(..., along=4),
        .packages=c("pomp", "spatPomp"),
        .options.multicore=mcopts) %dopar%
        {
          mult_rep_output[[i]]@prev_weights
        }
      prev_weights <- prev_prev_weights[,(dim(prev_prev_weights)[2]-max_lookback+1):dim(prev_prev_weights)[2],,,drop=FALSE]
      #  for log-likelihood computation
      cond_loglik[,nt] <- foreach::foreach(u=seq_len(nunits),
                                    .combine = 'c',
                                    .packages=c("pomp", "spatPomp"),
                                    .options.multicore=mcopts) %dopar%
      {
        log_mp_sum = logmeanexp(vapply(mult_rep_output,
                                       FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,]),
                                       FUN.VALUE = 1.0))
        log_p_sum = logmeanexp(vapply(mult_rep_output,
                                      FUN = function(rep_output) return(rep_output@log_wp_avg[u,]),
                                      FUN.VALUE = 1.0))
        cond_loglik_u = log_mp_sum - log_p_sum
        cond_loglik_u
      }
      #  for parameter swarm
      param_swarm <- foreach::foreach(
        i=seq_len(Nrep),
        .combine = 'cbind',
        .packages=c("pomp", "spatPomp"),
        .options.multicore=mcopts) %dopar%
        {
          mult_rep_output[[i]]@param
        }
      # parameter selection
      rep_param_init <- apply(param_swarm,1,mean)

      #  states
      states <- foreach::foreach(
        i=seq_len(Nrep),
        .combine = 'cbind',
        .packages=c("pomp", "spatPomp"),
        .options.multicore=mcopts) %dopar%
        {
          mult_rep_output[[i]]@state
        }

      rm(param_swarm)
      rm(mult_rep_output)
    } # end time loop
    # end single threaded

    # begin multi-threaded
    # pmag_init <- cooling.fn(1,n)$alpha*rw.sd[,1]*2
    # rep_param_init <- .Call('randwalk_perturbation',rep_param_init,pmag_init,PACKAGE = 'pomp')
    # param_swarm <- rep_param_init*0
    # num_batches <- ceiling(Nrep/reps_per_batch)
    # reps_in_batch <- c(rep(reps_per_batch, Nrep%/%reps_per_batch),Nrep%%reps_per_batch)
    # if(Nrep%%reps_per_batch != 0) {
    #   cum_reps_in_batch <- c(0,cumsum(reps_in_batch))
    #   batch_weight <- reps_in_batch/Nrep
    # } else {
    #   # no need for last element of reps_per_batch since it is 0
    #   cum_reps_in_batch <- c(0,cumsum(reps_in_batch)[-length(reps_in_batch)])
    #   batch_weight <- reps_in_batch[-length(reps_in_batch)]/Nrep
    # }
    # cond_loglik_un_batch <- array(dim=c(nunits,ntimes,num_batches))
    # loglik_rep <- vector(length=Nrep)
    # for(k in seq_len(num_batches)){
    #   mult_rep_output <- foreach::foreach(i=(cum_reps_in_batch[k]+1):(cum_reps_in_batch[k+1]),
    #                                       .packages=c("pomp","spatPomp"),
    #                                       .options.multicore=mcopts) %dopar%  {
    #                                         spatPomp:::h_abf_internal(
    #                                           object=object,
    #                                           params=rep_param_init[,i],
    #                                           Np=Np,
    #                                           nbhd=nbhd,
    #                                           abfiter=.ndone+n,
    #                                           cooling.fn=cooling.fn,
    #                                           rw.sd=rw.sd,
    #                                           tol=tol,
    #                                           max.fail=max.fail,
    #                                           .indices=.indices,
    #                                           verbose=verbose
    #                                         )
    #                                       }
    #   ## for log-likelihood computation
    #   cond_loglik <- foreach::foreach(u=seq_len(nunits),
    #                                   .combine = 'rbind',
    #                                   .packages=c("pomp", "spatPomp"),
    #                                   .options.multicore=mcopts) %dopar%
    #     {
    #       cond_loglik_u <- array(data = numeric(0), dim=c(ntimes))
    #       for (n in seq_len(ntimes)){
    #         log_mp_sum = logmeanexp(vapply(mult_rep_output,
    #                                        FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,n]),
    #                                        FUN.VALUE = 1.0))
    #         log_p_sum = logmeanexp(vapply(mult_rep_output,
    #                                       FUN = function(rep_output) return(rep_output@log_wp_avg[u,n]),
    #                                       FUN.VALUE = 1.0))
    #         cond_loglik_u[n] = log_mp_sum - log_p_sum
    #       }
    #       cond_loglik_u
    #     }
    #   cond_loglik_un_batch[,,k] <- cond_loglik
    #   rm(cond_loglik)
    #
    #   # for parameter estimation
    #   rep_loglik_un <- foreach::foreach(u=seq_len(nunits),
    #                                     .combine = function(...) abind::abind(..., along=3),
    #                                     .packages=c("pomp", "spatPomp"),
    #                                     .options.multicore=mcopts) %dopar%
    #     {
    #       cond_loglik_u <- array(data = numeric(0), dim=c(ntimes,length((cum_reps_in_batch[k]+1):(cum_reps_in_batch[k+1]))))
    #       for (n in seq_len(ntimes)){
    #         rep_filt_weight_un_rep = vapply(mult_rep_output,
    #                                         FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,n]),
    #                                         FUN.VALUE = 1.0) -
    #           vapply(mult_rep_output,
    #                  FUN = function(rep_output) return(rep_output@log_wp_avg[u,n]),
    #                  FUN.VALUE = 1.0)
    #         cond_loglik_u[n,] = rep_filt_weight_un_rep
    #       }
    #       cond_loglik_u
    #     }
    #   loglik_rep[(cum_reps_in_batch[k]+1):(cum_reps_in_batch[k+1])] <- apply(rep_loglik_un, MARGIN = 2, FUN = sum)
    #   rm(rep_loglik_un)
    #
    #   # for parameter swarm
    #   param_swarm[,(cum_reps_in_batch[k]+1):(cum_reps_in_batch[k+1])] <- foreach::foreach(
    #     i=seq_len(length((cum_reps_in_batch[k]+1):(cum_reps_in_batch[k+1]))),
    #     .combine = 'cbind',
    #     .packages=c("pomp", "spatPomp"),
    #     .options.multicore=mcopts) %dopar%
    #     {
    #       mult_rep_output[[i]]@param
    #     }
    # }
    # end multi-threaded
    gnsi <- FALSE
    coef(object) <- partrans(object,
                             rep_param_init,
                             dir="fromEst", .gnsi=.gnsi)

    traces[n+1,-c(1)] <- coef(object)
    traces[n+1,c(1)] <- sum(cond_loglik)
    if (verbose) {
      cat("iabf iteration",n,"of",Nabf,"completed","with log-likelihood",sum(cond_loglik),"\n")
      print(coef(object))
    }
  }

  # parameter swarm to be outputted
  param_swarm <- matrix(coef(object),nrow=length(coef(object)),ncol=1,dimnames=list(params=names(coef(object))))
  pompUnload(object,verbose=FALSE)

  new(
    "iabf2d_spatPomp",
    object,
    Nabf=as.integer(Nabf),
    Nrep=as.integer(Nrep),
    Np=as.integer(Np),
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces,
    paramMatrix=param_swarm,
    loglik=sum(cond_loglik)
  )
}

setGeneric(
  "iabf2",
  function (object, ...)
    standardGeneric("iabf2")
)

##' @name iabf2-iabf2d_spatPomp
##' @aliases iabf2,iabf2d_spatPomp-method
##' @rdname iabf2
##' @export
setMethod(
  "iabf2",
  signature=signature(object="spatPomp"),
  definition = function (object, Nabf = 1, Nrep, nbhd, Np,
                         rw.sd,
                         cooling.type = c("geometric","hyperbolic"),
                         cooling.fraction.50, tol = (1e-18)^17,
                         max.fail = Inf,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("iabf2"),": ")
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
    if (missing(rw.sd))
      stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
    cooling.type <- match.arg(cooling.type)
    cooling.fraction.50 <- as.numeric(cooling.fraction.50)
    if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
      stop(ep,sQuote("cooling.fraction.50"),
           " must be in (0,1]",call.=FALSE)

    iabf_internal2(
      object=object,
      Nrep=Nrep,
      Nabf=Nabf,
      nbhd=nbhd,
      Np=Np,
      rw.sd=rw.sd,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      ...
    )
  }
)
