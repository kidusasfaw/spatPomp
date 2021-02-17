##' Iterated Adapted Bagged Filter 6 (IABF6)
##'
##' An algorithm for estimating the parameters of a spatiotemporal partially-observed Markov process.
##' Running \code{iabf6} causes the algorithm to perform a specified number of iterations of adapted simulations with parameter perturbation and parameter resamplings.
##' At each iteration, adapted simulations are performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' After cycling through the data, each replicate's weight is calculated and is used to rank the bootstrap replictates. The highest ranking replicates are recycled into the next iteration.
##' This extra variability introduced through parameter perturbation effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##'
##' @name iabf6
##' @rdname iabf6
##' @include spatPomp_class.R abf.R
##' @family particle filter methods
##' @family \pkg{spatPomp} parameter estimation methods
##' @importFrom stats quantile
##' @importFrom utils head
##' @inheritParams pomp::mif2
##' @inheritParams abf
##' @param Nabf The number of iterations to perform
##' @param Nparam The number of parameters that will undergo the iterated perturbation
##' @param Nrep_per_param The number of replicates used to estimate the likelihood at a parameter
##' @param prop A numeric between 0 and 1. The top \code{prop}*100\% of the parameters are resampled at each observation
##' @return
##' Upon successful completion, \code{iabf} returns an object of class
##' \sQuote{iabfd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ extracts the point estimate }
##' }
NULL

rw.sd <- pomp:::safecall

setClass(
  "iabf6_iter",
  slots=c(
    cond_loglik = 'numeric',
    paramMatrix = 'matrix'
  ),
  prototype=prototype(
    cond_loglik = numeric(0),
    paramMatrix = matrix(0)
  )
)

## define the iabfd_spatPomp class
setClass(
  'iabf6d_spatPomp',
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


iabf_abf <- function (object,
                             params,
                             Nrep_per_param,
                             Np,
                             nbhd,
                             prop,
                             abfiter,
                             cooling.fn,
                             rw.sd,
                             tol = (1e-18)^17,
                             .indices = integer(0), verbose,
                             .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("iabf_abf"),": ")
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  abfiter <- as.integer(abfiter)
  Np <- as.integer(Np)
  Nparam <- dim(params)[2]
  Nrep_per_param <- as.integer(Nrep_per_param)
  Nislands <- as.integer(Nparam*Nrep_per_param)

  all_times <- time(object,t0=TRUE)
  ntimes <- length(all_times)-1
  nunits <- length(unit_names(object))
  cond_loglik <- vector(length=ntimes)
  prev_meas_weights <- NULL
  pompLoad(object,verbose=FALSE)
  resample_ixs_raw <- rep(1:Nparam)

  for(nt in seq_len(ntimes)){
    if(verbose && nt %in% c(1,2)) {
      cat("working on observation times ", nt, " in iteration ", abfiter, "at time ", Sys.time(), "\n")
    }
    if(verbose){
      print(c(min(pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)['psi',]),
              max(pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)['psi',])))
    }
    params <- params[,resample_ixs_raw]
    pmag <- cooling.fn(nt,abfiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')

    all_params <- params[,rep(1:Nparam, each = Np[1L]*Nrep_per_param)]
    tparams <- pomp::partrans(object,all_params,dir="fromEst",.gnsi=gnsi)

    gc()

    if (nt == 1L) X <- rinit(object,params=tparams)
    else X <- X[,resample_ixs]

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        x0=X,
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

    gc()

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
    log_cond_densities <- log_weights
    log_resamp_weights <- colSums(log_weights[,,1])
    max_positions <- max.col(t(matrix(log_resamp_weights, ncol = Nislands)))
    max_positions_vec <- ((1:Nislands)-1)*Np[1L] + max_positions
    max_positions_vec <- rep(max_positions_vec, each = Np[1L])
    log_resamp_weights <- log_resamp_weights - log_resamp_weights[max_positions_vec]
    # if any particle's resampling weight is zero replace by tolerance
    log_resamp_weights[is.infinite(log_resamp_weights)] <- log(tol)
    resamp_weights <- exp(log_resamp_weights)
    resamp_weights <- matrix(resamp_weights, nrow=Np[1L])
    u <- runif(ncol(resamp_weights))
    cumul_w <- lower.tri(diag(nrow(resamp_weights)), diag = TRUE)%*%(resamp_weights)
    cumul_w_norm <- cumul_w/colSums(resamp_weights)[col(cumul_w)]
    i <- colSums(u > cumul_w_norm) + 1L
    choices <- 1:Np[1L]
    selections <- choices[i]
    X <- X[,rep(selections + seq(from=0, to=((Nparam*Nrep_per_param)-1)*Np[1L], by = Np[1L]),each=Np[1L]),]
    rm(log_weights, log_resamp_weights, max_positions,max_positions_vec, resamp_weights, u, cumul_w, i, choices, selections)
    gnsi <- FALSE

    gc()

    if(!is.null(prev_meas_weights)) num_old_times <- dim(prev_meas_weights)[3]
    else num_old_times <- 0
    prev_meas_weights <- prev_meas_weights
    log_loc_comb_pred_weights = array(data = numeric(0), dim=c(nunits,Np[1L], Nrep_per_param*Nparam))
    log_wm_times_wp_avg = array(data = numeric(0), dim = c(nunits, Nrep_per_param*Nparam))
    log_wp_avg = array(data = numeric(0), dim = c(nunits, Nrep_per_param*Nparam))
    for (unit in seq_len(nunits)){
      full_nbhd <- nbhd(object, time = nt, unit = unit)
      log_prod_cond_dens_nt  <- rep(0, Np[1]*Nrep_per_param*Nparam)
      farthest_time <- nt-num_old_times
      log_prod_cond_dens_not_nt <- matrix(0, Np[1]*Nrep_per_param*Nparam, max(1,nt-farthest_time))
      for (neighbor in full_nbhd){
        neighbor_u <- neighbor[1]
        neighbor_n <- neighbor[2]
        if (neighbor_n == nt)
          log_prod_cond_dens_nt  <- log_prod_cond_dens_nt + log_cond_densities[neighbor_u, ,1]
        else{
          # means prev_meas_weights was non-null, i.e. dim(prev_meas_weights)[3]>=1
          log_prod_cond_dens_not_nt[, nt-neighbor_n] <- log_prod_cond_dens_not_nt[, nt-neighbor_n] +
            prev_meas_weights[neighbor_u, ,num_old_times+1-(nt-neighbor_n)]
        }
      }
      log_prod_cond_dens_not_nt_by_island <- array(log_prod_cond_dens_not_nt, dim = c(Np[1], Nrep_per_param*Nparam, max(1,num_old_times)))
      log_prod_cond_dens_nt_by_island <- matrix(log_prod_cond_dens_nt, nrow = Np[1L], ncol = Nrep_per_param*Nparam)
      first_term <- rowSums(apply(log_prod_cond_dens_not_nt_by_island, c(2,3), logmeanexp))
      second_term <- log_prod_cond_dens_nt_by_island
      log_loc_comb_pred_weights[unit,,]  <- first_term[col(second_term)] + second_term
    }

    log_cond_densities_by_island <- array(log_cond_densities, dim=c(nunits, Np[1L], Nrep_per_param*Nparam))
    log_wm_times_wp_avgs = apply(log_loc_comb_pred_weights + log_cond_densities_by_island, c(1,3), FUN = logmeanexp)
    log_wm_times_wp_avgs_by_param = apply(array(log_wm_times_wp_avgs, dim = c(nunits, Nrep_per_param, Nparam)), c(1,3), FUN = logmeanexp)

    log_wp_avgs = apply(log_loc_comb_pred_weights, c(1,3), FUN = logmeanexp)
    log_wp_avgs_by_param = apply(array(log_wp_avgs, dim = c(nunits, Nrep_per_param, Nparam)), c(1,3), FUN = logmeanexp)

    param_resamp_log_weights <- colSums(log_wm_times_wp_avgs_by_param - log_wp_avgs_by_param)

    ####### Quantile resampling
    def_resample <- which(param_resamp_log_weights > quantile(param_resamp_log_weights, 1-prop))
    length_also_resample <- Nparam - length(def_resample)
    also_resample <- sample(def_resample, size = length_also_resample, replace = TRUE, prob = rep(1, length(def_resample)))
    resample_ixs_raw <- c(def_resample, also_resample)
    rm(def_resample,also_resample, length_also_resample)
    ####### Resampling using parameter log-likelihood
    # param_resamp_weights <- exp(param_resamp_log_weights - max(param_resamp_log_weights))
    # resample_ixs_raw <- sample(1:Nparam, size = Nparam, replace = TRUE, prob = param_resamp_weights)
    resample_ixs <- (resample_ixs_raw - 1)*Np[1L]*Nrep_per_param
    resample_ixs <- rep(resample_ixs, each = Np[1L]*Nrep_per_param)
    resample_ixs <- rep(1:(Np[1L]*Nrep_per_param), Nparam) + resample_ixs

    # for next observation time, how far back do we need
    # to provide the conditional densities, f_{Y_{u,n}|X_{u,n}}?
    max_lookback <- 0
    for(u in seq_len(nunits)){
      all_nbhd_times <- sapply(nbhd(object=object,unit=u,time=nt+1),'[[',2)
      if(length(all_nbhd_times) == 0) next
      smallest_nbhd_time <- min(all_nbhd_times)
      if(nt+1-smallest_nbhd_time > max_lookback) max_lookback <- nt+1-smallest_nbhd_time
    }
    if(max_lookback == 1) prev_meas_weights <- log_cond_densities[,resample_ixs,1,drop=FALSE]
    if(max_lookback > 1){
      prev_meas_weights <- abind::abind(prev_meas_weights[,,(dim(prev_meas_weights)[3]+2-max_lookback):dim(prev_meas_weights)[3],drop=F],
                                        log_cond_densities,
                                        along=3)[,resample_ixs,,drop=FALSE]
    }
    cond_loglik[nt] <- logmeanexp(param_resamp_log_weights)
    gc()
  }
  params <- params[,resample_ixs_raw]
  pompUnload(object,verbose=FALSE)
  new(
    "iabf6_iter",
    cond_loglik = cond_loglik,
    paramMatrix=params
  )
}

iabf_internal6 <- function (object, Nrep_per_param, Nparam, nbhd, Nabf, Np, prop, rw.sd,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17,
                           verbose = FALSE, .ndone = 0L,
                           .indices = integer(0),
                           .paramMatrix = NULL,
                           .gnsi = TRUE, ...) {

  ep <- paste0("in ",sQuote("iabf6"),": ")
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
  rep_param_init <- array(data = rep_param_init,
                          dim = c(length(rep_param_init), Nparam),
                          dimnames = list(param = names(rep_param_init), rep = NULL))


  # iterate the filtering
  times = time(object)
  ntimes = length(time(object))
  nunits = length(unit_names(object))
  for (n in seq_len(Nabf)) {
    out <- iabf_abf(
      object=object,
      params=rep_param_init,
      Nrep_per_param=Nrep_per_param,
      Np=Np,
      nbhd=nbhd,
      prop=prop,
      abfiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      tol=tol,
      .indices=.indices,
      verbose=verbose
    )
    gnsi <- FALSE
    rep_param_init <- out@paramMatrix
    out_params_summary <- apply(partrans(object,
                                         rep_param_init,
                                         dir="fromEst", .gnsi=.gnsi),1,mean)

    traces[n+1,-c(1)] <- out_params_summary
    traces[n+1,c(1)] <- sum(out@cond_loglik)
    coef(object) <- out_params_summary
    if (verbose) {
      cat("iabf iteration",n,"of",Nabf,"completed","with log-likelihood",sum(out@cond_loglik),"\n")
      print(out_params_summary)
    }
  }
  # parameter swarm to be outputted
  param_swarm <- partrans(object,
                          rep_param_init,
                          dir="fromEst", .gnsi=.gnsi)
  pompUnload(object,verbose=FALSE)

  new(
    "iabf6d_spatPomp",
    object,
    Nabf=as.integer(Nabf),
    Nrep=as.integer(Nrep_per_param),
    Np=as.integer(Np),
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces,
    paramMatrix=param_swarm,
    loglik=sum(out@cond_loglik)
  )

}

setGeneric(
  "iabf6",
  function (object, ...)
    standardGeneric("iabf6")
)

##' @name iabf6-spatPomp
##' @aliases iabf6,spatPomp-method
##' @rdname iabf6
##' @export
setMethod(
  "iabf6",
  signature=signature(object="spatPomp"),
  definition = function (object, Nabf = 1, Nrep_per_param, Nparam, nbhd, Np,prop,
                         rw.sd,
                         cooling.type = c("geometric","hyperbolic"),
                         cooling.fraction.50, tol = (1e-18)^17,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("iabf6"),": ")
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

    iabf_internal6(
      object=object,
      Nrep_per_param=Nrep_per_param,
      Nparam=Nparam,
      Nabf=Nabf,
      nbhd=nbhd,
      Np=Np,
      prop=prop,
      rw.sd=rw.sd,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      verbose=verbose,
      ...
    )
  }
)
