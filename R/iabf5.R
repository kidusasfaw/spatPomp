##' Iterated Adapted Bagged Filter 5 (IABF5)
##'
##' An algorithm for estimating the parameters of a spatiotemporal partially-observed Markov process.
##' Running \code{iabf} causes the algorithm to perform a specified number of iterations of adapted simulations with parameter perturbation and parameter resamplings.
##' At each iteration, adapted simulations are performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' After cycling through the data, each replicate's weight is calculated and is used to rank the bootstrap replictates. The highest ranking replicates are recycled into the next iteration.
##' This extra variability introduced through parameter perturbation effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##'
##' @name iabf5
##' @rdname iabf5
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

setClass(
  "adapted_replicate_extended3",
  contains='adapted_replicate',
  slots=c(
    param = 'numeric',
    state = 'numeric',
    prev_weights = 'array',
    island_weight = "numeric"
  ),
  prototype=prototype(
    param=numeric(0),
    state=numeric(0),
    prev_weights = array(data=numeric(0),dim=c(0,0,0)),
    island_weight = numeric(0)
  )
)


h_abf_internal5 <- function (object,
                             params,
                             states,
                             obs_nums,
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

  all_times <- time(object,t0=TRUE)
  n_obs_times <- length(obs_nums)
  ntimes <- length(all_times)-1
  nunits <- length(unit_names(object))

  params <- array(data = params,
                  dim = c(length(params),Np[1L]),
                  dimnames = list(param = names(params), rep = NULL))

  if(!is.null(prev_meas_weights)){
    prev_meas_weights <- array(data = prev_meas_weights,
                               dim = c(dim(prev_meas_weights)[1],
                                       dim(prev_meas_weights)[2],
                                       dim(prev_meas_weights)[3]),
                               )
  }
  if(!is.null(states)) x <- states
  if(!is.null(prev_meas_weights)) num_old_times <- dim(prev_meas_weights)[3]
  else num_old_times <- 0
  num_new_times <- 0
  # create array to store weights across time
  log_cond_densities <- array(data = numeric(0), dim=c(nunits,Np[1L],n_obs_times))
  dimnames(log_cond_densities) <- list(unit = 1:nunits, rep = 1:Np[1L], time = obs_nums)
  log_cond_densities <- abind::abind(prev_meas_weights, log_cond_densities, along=3)
  rm(prev_meas_weights)
  pompLoad(object,verbose=FALSE)
  for(nt in obs_nums){
    num_new_times <- num_new_times+1
    # NO PARAMETER PERTURBATIONS INSIDE OF A REPLICATE
    tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)
    if (nt == 1L) x <- rinit(object,params=tparams)

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        x0=x,
        t0=all_times[nt],
        times=all_times[nt+1],
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
        times=all_times[nt+1],
        params=tparams,
        log=TRUE,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    log_cond_densities[,,num_old_times+num_new_times] <- log_weights[,,1]
    log_resamp_weights <- apply(log_weights[,,1,drop=FALSE], 2, function(x) sum(x))
    max_log_resamp_weights <- max(log_resamp_weights)
    # if any particle's resampling weight is zero replace by tolerance
    if(all(is.infinite(log_resamp_weights))) log_resamp_weights <- rep(log(tol), Np[1L])
    else log_resamp_weights <- log_resamp_weights - max_log_resamp_weights
    resamp_weights <- exp(log_resamp_weights)
    gnsi <- FALSE
    xx <- tryCatch(
      .Call(
        "abf_computations",
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
  }

  log_loc_comb_pred_weights = array(data = numeric(0), dim=c(nunits,Np[1L], n_obs_times))
  log_wm_times_wp_avg = array(data = numeric(0), dim = c(nunits, n_obs_times))
  log_wp_avg = array(data = numeric(0), dim = c(nunits, n_obs_times))
  num_new_times <- 0
  for (nt in obs_nums){
    num_new_times <- num_new_times + 1
    for (unit in seq_len(nunits)){
      full_nbhd <- nbhd(object, time = nt, unit = unit)
      log_prod_cond_dens_nt  <- rep(0, Np[1])
      if(length(full_nbhd) > 0) log_prod_cond_dens_not_nt <- matrix(0, Np[1], max(1,nt-min(sapply(full_nbhd,'[[',2))))
      else log_prod_cond_dens_not_nt <- matrix(0,Np[1],0)
      for (neighbor in full_nbhd){
        neighbor_u <- neighbor[1]
        neighbor_n <- neighbor[2]
        if (neighbor_n == nt)
          log_prod_cond_dens_nt  <- log_prod_cond_dens_nt + log_cond_densities[neighbor_u, ,num_old_times + num_new_times]
        else{
          # means prev_meas_weights was non-null
          if(dim(log_cond_densities)[3]>length(obs_nums)){
            log_prod_cond_dens_not_nt[, nt-neighbor_n] <- log_prod_cond_dens_not_nt[, nt-neighbor_n] +
              log_cond_densities[neighbor_u,,num_old_times+num_new_times-(nt-neighbor_n)]
          } else {
            log_prod_cond_dens_not_nt[, nt-neighbor_n] <- log_prod_cond_dens_not_nt[, nt-neighbor_n] +
              log_cond_densities[neighbor_u, ,neighbor_n]
          }
        }
      }
      log_loc_comb_pred_weights[unit, ,num_new_times]  <- sum(apply(log_prod_cond_dens_not_nt, 2, logmeanexp)) + log_prod_cond_dens_nt
    }
  }
  # if(dim(log_cond_densities)[3]>length(obs_nums))
  #   # log_cond_densities comes as a U*J*N so flip it back to U*N*J for iabf_internal
  #   prev_meas_weights <- array(log_cond_densities, dim=c(dim(log_cond_densities)[1],
  #                                                        dim(log_cond_densities)[3],
  #                                                        dim(log_cond_densities)[2]))
  # else
  #   prev_meas_weights <- array(log_cond_densities, dim=c(nunits, num_new_times, Np[1]))
  prev_meas_weights <- log_cond_densities
  params_last = params[,1]
  log_wm_times_wp_avg = apply(log_loc_comb_pred_weights + log_cond_densities[,,(num_old_times+1):
                                                                               (dim(log_cond_densities)[3]),drop=FALSE], c(1,3), FUN = logmeanexp)
  log_wp_avg = apply(log_loc_comb_pred_weights, c(1,3), FUN = logmeanexp)
  pompUnload(object,verbose=FALSE)
  new(
    "adapted_replicate_extended3",
    log_wm_times_wp_avg = log_wm_times_wp_avg,
    log_wp_avg = log_wp_avg,
    param=params_last,
    state=x[,1],
    prev_weights=prev_meas_weights,
    island_weight = sum(log_wm_times_wp_avg - log_wp_avg),
    Np=as.integer(Np),
    tol=tol
  )
}

iabf_internal5 <- function (object, Nrep_per_param, Nparam, nbhd, Nabf, Np, resample_every, prop, parallel, rw.sd,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17, max.fail = Inf,
                           verbose = FALSE, .ndone = 0L,
                           .indices = integer(0),
                           .paramMatrix = NULL,
                           .gnsi = TRUE, ...) {

  ep <- paste0("in ",sQuote("iabf5"),": ")
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
  resample_every <- as.integer(resample_every)

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
  rep_param_init <- array(data = rep_param_init,
                          dim = c(length(rep_param_init), Nparam),
                          dimnames = list(param = names(rep_param_init), rep = NULL))


  # iterate the filtering
  mcopts <- list(set.seed=TRUE)
  times = time(object)
  ntimes = length(time(object))
  nunits = length(unit_names(object))
  obs_seq = seq_len(ntimes)
  nblocks = round(length(obs_seq)/resample_every)
  block_list = split(obs_seq, sort(obs_seq %% nblocks))
  for (n in seq_len(Nabf)) {
    prev_weights <- NULL
    cond_loglik <- array(dim=c(nunits,ntimes))
    states <- NULL
    for (b in block_list){
      if(verbose) cat("working on observation times ", b, " in iteration ", n, "\n")
      mult_rep_output <- list()
      ## perturb parameters
      pmag <- cooling.fn(min(b),n)$alpha*rw.sd[,min(b)]
      rep_param_init <- .Call('randwalk_perturbation',rep_param_init,pmag,PACKAGE = 'pomp')
      if(Nparam > 50) print(mean(partrans(object,rep_param_init,dir="fromEst",.gnsi=gnsi)['psi',]))
      all_params <- rep_param_init[,rep(1:Nparam, each=Nrep_per_param)]
      if(parallel){
        # begin multi-threaded
        mult_rep_output <- foreach::foreach(i=seq_len(Nparam*Nrep_per_param),
                                            .packages=c("pomp","spatPomp"),
                                            .options.multicore=mcopts) %dopar%  {
                                              spatPomp:::h_abf_internal5(
                                                object=object,
                                                params=all_params[,i],
                                                states=(if(is.null(states)) NULL else states[,i]),
                                                obs_nums=b,
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
                                              )
                                            }
        # end multi-threaded
      } else{
        # begin single-threaded
        for(i in seq_len(Nparam*Nrep_per_param)){
          mult_rep_output <- c(mult_rep_output, spatPomp:::h_abf_internal5(
            object=object,
            params=all_params[,i],
            states=(if(is.null(states)) NULL else states[,i]),
            obs_nums=b,
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
        # end single-threaded
      }
      # for the next observation time, how far back do we need
      # to provide the conditional densities, f_{Y_{u,n}|X_{u,n}}?
      max_lookback <- 0
      for(u in seq_len(nunits)){
        farthest_back <- min(sapply(nbhd(unit=u,time=(max(b)+1)),'[[',2))
        if(max(b)+1 - farthest_back > max_lookback) max_lookback <- max(b)+1 - farthest_back
      }
      # update prev_weights for the next observation time
      # THIS NEEDS TIGHTENING UP FOR DIFFERENT KINDS OF NEIGHBORHOODS
      prev_prev_weights <- foreach::foreach(
        i=seq_len(Nparam*Nrep_per_param),
        .combine = function(...) abind::abind(..., along=4),
        .packages=c("pomp", "spatPomp"),
        .options.multicore=mcopts) %dopar%
        {
          mult_rep_output[[i]]@prev_weights
        }
      prev_weights <- prev_prev_weights[,,(dim(prev_prev_weights)[3]-max_lookback+1):dim(prev_prev_weights)[3],,drop=FALSE]
      #  for log-likelihood computation
      cond_loglik[,b] <- foreach::foreach(u=seq_len(nunits),
                                      .combine = 'rbind',
                                      .packages=c("pomp", "spatPomp"),
                                      .options.multicore=mcopts) %dopar%
        {
          cond_loglik_u <- array(data = numeric(0), dim=c(length(b)))
          for (n in seq_len(length(b))){
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
      #  for param_weights. think i'm calculating parameter weights for all times in b and then adding them up???
      param_weights <- foreach::foreach(u=seq_len(nunits),
                                          .combine = function(...) abind::abind(..., along=3),
                                          .packages=c("pomp", "spatPomp"),
                                          .options.multicore=mcopts) %dopar%
        {
          param_weight_u <- array(data = numeric(0), dim=c(length(b), Nparam))
          for (n in seq_len(length(b))){
            all_log_wm_times_wp_avg <- vapply(mult_rep_output,
                                              FUN = function(rep_output) return(rep_output@log_wm_times_wp_avg[u,n]),
                                              FUN.VALUE = 1.0)
            all_log_wp_avg <- vapply(mult_rep_output,
                                 FUN = function(rep_output) return(rep_output@log_wp_avg[u,n]),
                                 FUN.VALUE = 1.0)

            all_log_wm_times_wp_avg_matrix <- matrix(all_log_wm_times_wp_avg, nrow = Nrep_per_param)
            all_log_wp_avg_matrix <- matrix(all_log_wp_avg, nrow = Nrep_per_param)
            avg_log_wm_times_wp_avg_by_param <- apply(all_log_wm_times_wp_avg_matrix, 2, logmeanexp)
            avg_log_wp_avg_by_param <- apply(all_log_wp_avg_matrix, 2, logmeanexp)
            param_weight_u[n,] = avg_log_wm_times_wp_avg_by_param - avg_log_wp_avg_by_param
          }
          param_weight_u
        }
      loglik_param <- apply(param_weights, 2, FUN = sum)
      # parameter selection
      def_resample <- which(loglik_param > quantile(loglik_param, 1-prop))
      length_also_resample <- Nparam - length(def_resample)
      also_resample <- sample(def_resample, size = length_also_resample, replace = TRUE, prob = rep(1, length(def_resample)))
      resampled_idx <- c(def_resample, also_resample)
      rep_param_init <- rep_param_init[,resampled_idx]
      rm(def_resample,also_resample, length_also_resample)

      #  states
      states <- lapply(mult_rep_output, slot, name="state")
      states <- do.call(cbind, states)
      resampled_idx2 <- rep(resampled_idx, each = Nrep_per_param)
      resampled_idx2_mat <- matrix(resampled_idx2, ncol = Nparam)
      resampled_idx3_mat <- (resampled_idx2_mat - 1)*Nrep_per_param
      to_add_mat <- matrix(rep(1:Nrep_per_param, Nparam), nrow = Nrep_per_param)
      resampled_idx_final <- as.integer(resampled_idx3_mat + to_add_mat)
      states <- states[,resampled_idx_final]
      prev_weights <- prev_weights[,,,resampled_idx_final,drop=FALSE]
      rm(resampled_idx, resampled_idx2, resampled_idx2_mat, resampled_idx3_mat, to_add_mat, resampled_idx_final)
      rm(mult_rep_output)
    } # end block loop
    gnsi <- FALSE
    coef(object) <- apply(partrans(object,
                             rep_param_init,
                             dir="fromEst", .gnsi=.gnsi),1,mean)

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
    Nrep=as.integer(Nrep_per_param),
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
  "iabf5",
  function (object, ...)
    standardGeneric("iabf5")
)

##' @name iabf5
##' @aliases iabf5
##' @rdname iabf5
##' @export
setMethod(
  "iabf5",
  signature=signature(object="spatPomp"),
  definition = function (object, Nabf = 1, Nrep_per_param, Nparam, nbhd, Np,resample_every,prop,
                         parallel = TRUE, rw.sd,
                         cooling.type = c("geometric","hyperbolic"),
                         cooling.fraction.50, tol = (1e-18)^17,
                         max.fail = Inf,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("iabf5"),": ")
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
    if(missing(Nrep_per_param))
      stop(ep,"number of replicates, ",
           sQuote("Nrep_per_param"),", must be specified!",call.=FALSE)
    if (missing(rw.sd))
      stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
    cooling.type <- match.arg(cooling.type)
    cooling.fraction.50 <- as.numeric(cooling.fraction.50)
    if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
      stop(ep,sQuote("cooling.fraction.50"),
           " must be in (0,1]",call.=FALSE)

    iabf_internal5(
      object=object,
      Nrep_per_param=Nrep_per_param,
      Nparam=Nparam,
      Nabf=Nabf,
      nbhd=nbhd,
      Np=Np,
      resample_every=resample_every,
      prop=prop,
      parallel=parallel,
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
