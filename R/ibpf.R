##' Iterated block particle filter (IBPF)
##'
##' An iterated block particle filter, for both shared and unit-specific
##' parameters. We require that the spatPomp has
##' been constructed to have a unit-specific parameter "thetau"
##' for unit u corresponding to an estimated parameter "theta", whether
##' theta is shared or unit-specific. This permits IBPF
##' to implement a spatiotemporal random walk to estimate theta.
##' We require that rw.sd is positive for, and only for, all parameters
##' of the form "thetau" if "theta" is listed in sharedParNames or
##' unitParNames.
##'
##' @name ibpf
##' @rdname ibpf
##' @include spatPomp_class.R spatPomp.R bpfilter.R iter_filter.R
##' @author Edward L. Ionides
##' @family likelihood maximization algorithms
##' @seealso likelihood evaluation algorithms: \code{girf()}, \code{enkf()}, \code{bpfilter()}, \code{abf()}, \code{abfir()}
##'
##' @importFrom foreach %do%
##'
##' @inheritParams bpfilter
##' @inheritParams pomp::mif2
##' @param Nbpf the number of iterations of perturbed BPF.
##' @param sharedParNames estimated parameters that are equal for each unit.
##' @param unitParNames estimated parameters that are different for
##' each unit.
##' @param spat_regression fraction of each extended parameter regressed toward the unit mean. Not required when all parameters are unit-specific.
##'
##' @return
##' Upon successful completion, \code{ibpf} returns an object of class
##' \sQuote{ibpfd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ gives the Monte Carlo estimate of the maximum likelihood. }
##' }
##'
##' @references \ionides2022
##'

NULL

rw.sd <- safecall

setClass(
  "ibpfd_spatPomp",
  contains="bpfilterd_spatPomp",
  slots=c(
    Nbpf = 'integer',
    spat_regression = 'numeric',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix'
  )
)

setGeneric(
  "ibpf",
  function (data, ...)
    standardGeneric("ibpf")
)

##' @name ibpf-missing
##' @aliases ibpf,missing-method
##' @rdname ibpf
##' @export
setMethod(
  "ibpf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("ibpf","data")
  }
)

##' @name ibpf-ANY
##' @aliases ibpf,ANY-method
##' @rdname ibpf
##' @export
setMethod(
  "ibpf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("ibpf",data)
  }
)

##' @name ibpf-spatPomp
##' @aliases ibpf,spatPomp-method
##' @rdname ibpf
##' @export
setMethod(
  "ibpf",
  signature=signature(data="spatPomp"),
  definition=function (data,Nbpf,Np,rw.sd,
    sharedParNames,
    unitParNames,
    cooling.type="geometric",
    cooling.fraction.50,
    block_size, block_list,spat_regression,
    ..., verbose = getOption("verbose", FALSE)
  ){
    ep <- paste0("in ",sQuote("ibpf"),": ")
    if (missing(Nbpf)) pStop_(ep, "Nbpf is required")
    if (missing(rw.sd)) pStop_(ep, "rw.sd is required")  
    if (missing(Np)) pStop_(ep, "Np is required")
    if (missing(sharedParNames)) pStop_(ep, "sharedParNames is required")
    if (missing(unitParNames)) pStop_(ep, "unitParNames is required")
    if (missing(spat_regression)) spat_regression <- numeric()
    if(missing(block_list) && missing(block_size)) {
      pStop_(ep,sQuote("block_list"), " or ", sQuote("block_size"),
        " must be specified to the call")
    }
    if (!missing(block_list) & !missing(block_size)){
      pStop_(ep,"Exactly one of ",sQuote("block_size"), " and ",
        sQuote("block_list"), " should be provided, but not both.")
    }
    if (missing(spat_regression) && !is.null(sharedParNames)) 
      pStop_(ep, sQuote("spat_regression"),
       " should be provided when there are shared parameters")
      if (missing(block_list)){
      if(block_size > length(unit_names(data))){
        pStop_(ep,sQuote("block_size"), " cannot be greater than the number of spatial units")
      }
      all_units = seq_len(length(unit_names(data)))
      nblocks = round(length(all_units)/block_size)
      block_list = split(all_units, sort(all_units %% nblocks))
    }
    block_list <- lapply(block_list, as.integer)
    if (missing(cooling.fraction.50)) pStop_(ep, "cooling.fraction.50 is required")

    tryCatch(
      ibpf_internal(data,
        Nbpf=Nbpf,
        spat_regression=spat_regression,
        Np=Np,rw.sd=rw.sd,
        sharedParNames=sharedParNames,
        unitParNames=unitParNames,
        cooling.type=cooling.type,
	cooling.fraction.50=cooling.fraction.50,
        block_list=block_list,
        ...,verbose=verbose
      ),
      error = function(e) pStop("ibpf",conditionMessage(e))
    )
  }
)


##' @name ibpf-ibpfd_spatPomp
##' @aliases ibpf,ibpfd_spatPomp-method
##' @rdname ibpf
##' @export
setMethod(
  "ibpf",
  signature=signature(data="ibpfd_spatPomp"),
  definition=function (data,Nbpf,Np,rw.sd,
    sharedParNames,
    unitParNames,
    cooling.type="geometric",
    cooling.fraction.50,
    block_size, block_list,spat_regression,
    ..., verbose = getOption("verbose", FALSE)
  ){
    ep <- paste0("in ",sQuote("ibpf"),": ")
    if (!missing(block_list) & !missing(block_size)){
      stop(ep,"Exactly one of ",sQuote("block_size"), " and ", sQuote("block_list"), " can be provided, but not both.",call.=FALSE)
    }

    if(missing(block_list) && missing(block_size))
      block_list <- data@block_list

    if (!missing(block_size)){
      if(block_size > length(unit_names(data))){
        stop(ep,sQuote("block_size"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      all_units = seq_len(length(unit_names(data)))
      nblocks = round(length(all_units)/block_size)
      block_list = split(all_units, sort(all_units %% nblocks))
    }
    block_list <- lapply(block_list, as.integer)

    if (missing(Np)) Np <- data@Np
    if (missing(Nbpf)) Nbpf <- data@Nbpf
    if (!is.numeric(Nbpf)|| Nbpf <1) stop(ep, sQuote("Nbpf"), " should be a positive integer")
    if (missing(rw.sd)) rw.sd <- data@rw.sd
    if (missing(cooling.type)) cooling.type <- data@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- data@cooling.fraction.50

    tryCatch(
      ibpf_internal(data,Nbpf=Nbpf,spat_regression=spat_regression,
        Np=Np,rw.sd=rw.sd,
	sharedParNames=sharedParNames,
	unitParNames=unitParNames,
	cooling.type=cooling.type,
	cooling.fraction.50=cooling.fraction.50,
        block_list=block_list,
	...,verbose=verbose),
        error = function(e) pStop("ibpf",conditionMessage(e))
    )
  }
)

##' @name ibpf-bpfd_spatPomp
##' @aliases ibpf,bpfilterd_spatPomp-method
##' @rdname ibpf
##' @export
setMethod(
  "ibpf",
  signature=signature(data="bpfilterd_spatPomp"),
  definition=function (data,Nbpf,Np,rw.sd,
    sharedParNames,
    unitParNames,
    cooling.type="geometric",
    cooling.fraction.50,
    block_size, block_list,spat_regression,
    ..., verbose = getOption("verbose", FALSE)
  ) {
    ep <- paste0("in ",sQuote("ibpf"),": ")
    if (!missing(block_list) & !missing(block_size)){
      stop(ep,"Exactly one of ",sQuote("block_size"), " and ", sQuote("block_list"), " can be provided, but not both.",call.=FALSE)
    }
    if(missing(block_list) && missing(block_size)) block_list <- data@block_list
    if (!missing(block_size)){
      if(block_size > length(unit_names(data))){
        stop(ep,sQuote("block_size"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      all_units = seq_len(length(unit_names(data)))
      nblocks = round(length(all_units)/block_size)
      block_list = split(all_units, sort(all_units %% nblocks))
    }
    block_list <- lapply(block_list, as.integer)

    if (missing(Np)) Np <- data@Np

    tryCatch(
      ibpf_internal(data,Nbpf=Nbpf,spat_regression=spat_regression,
        Np=Np,rw.sd=rw.sd,
	sharedParNames=sharedParNames,
	unitParNames=unitParNames,
	cooling.type=cooling.type,
	cooling.fraction.50=cooling.fraction.50,
        block_list=block_list,
	.paramMatrix=NULL,
	...,verbose=verbose),
      error = function (e) pStop(who="ibpf",conditionMessage(e))
    )
  }
)


ibpf_internal <- function (object,Nbpf,spat_regression,
   Np, rw.sd,cooling.type,cooling.fraction.50,
   sharedParNames,
   unitParNames,
   block_list, ...,
   .ndone = 0L, .indices = integer(0),.paramMatrix = NULL,
   .gnsi = TRUE, verbose = FALSE) {

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
  if (undefined(object@rprocess) || undefined(object@dunit_measure))
    pStop_(paste(sQuote(c("rprocess","dunit_measure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Nbpf) != 1 || !is.numeric(Nbpf) || !is.finite(Nbpf) || Nbpf < 1)
    pStop_(sQuote("Nbpf")," must be a positive integer.")
  Nbpf <- as.integer(Nbpf)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {  ## if '.paramMatrix' is supplied, 'start' is ignored
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))
  U <- length(unit_names(object))

  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
               " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1) {
    Np <- rep(Np,times=ntimes)
  } else if (length(Np) > ntimes) {
    if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
      pWarn("ibpf","Np[k] ignored for k > ",sQuote("length(time(object))"),".")
    }
    Np <- head(Np,ntimes)
  } else if (length(Np) < ntimes) {
    pStop_(sQuote("Np")," must have length 1 or ",
           sQuote("length(time(object))"),".")
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")

  Np <- as.integer(Np)
  Np <- c(Np,Np[1L])

  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
  
  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }

  traces <- array(dim=c(Nbpf+1,length(start)+1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Nbpf),
                                variable=c('loglik',names(start))))
  traces[1L,] <- c(NA,start)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)


  ## set up for the autoregression
  sharedParNamesExpanded <- sapply(sharedParNames,function(z) paste0(z,1:U))
  attr(sharedParNamesExpanded,"dim") <- NULL
  estParNames <- c(sharedParNames,unitParNames)
  
  ## iterate the filtering
  for (m in seq_len(Nbpf)) {
    b <- ibpf_bpfilter(object=object,
      block_list=block_list,params=paramMatrix,spat_regression=spat_regression,
      sharedParNames=sharedParNames,
      sharedParNamesExpanded=sharedParNamesExpanded,
      estParNames=estParNames,
      Np=Np,nbpf=.ndone+m,cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      verbose=verbose,.indices=.indices, .gnsi=gnsi)
    gnsi <- FALSE
    paramMatrix <- b@paramMatrix
    traces[m+1,-1] <- coef(b)
    traces[m+1,1] <- b@loglik
    .indices <- .indices

    if (verbose) cat("ibpf iteration",m,"of",Nbpf,"completed\n")

  }

  coef(b) <- mean_by_unit(coef(b),expandedParNames=sharedParNames,U=U)

  b@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
                            .gnsi=gnsi)

  new(
    "ibpfd_spatPomp",
    b,
    Nbpf=Nbpf,
    spat_regression=spat_regression,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

ibpf_bpfilter <- function (object,block_list,params,
    spat_regression,
    sharedParNames,
    sharedParNamesExpanded,
    estParNames,
    Np,nbpf,cooling.fn,rw.sd, 
    verbose,.indices = integer(0), .gnsi = TRUE) {

  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  nbpf <- as.integer(nbpf)
  Np <- as.integer(Np)

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  U <- length(unit_names(object))
  nblocks <- length(block_list)
  loglik <- rep(NA,ntimes)

  if (length(Np)==1) {
    Np <- rep(Np,times=ntimes+1)
  } else if (length(Np)!=(ntimes+1)) {
    pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1)
  } else if (!all(Np==Np[1]))
    pStop_("time-varying ", sQuote("Np")," is currently unsupported")
  if (any(Np<=0))
    pStop_("number of particles, ",sQuote("Np"),", must always be positive")
  if (!is.numeric(Np))
    pStop_(sQuote("Np"),
      " must be a number, a vector of numbers, or a function")
  Np <- as.integer(Np)
  if (is.matrix(params)) {
    if (!all(Np==ncol(params)))
      pStop_("when ",sQuote("params"),
        " is provided as a matrix, do not specify ", sQuote("Np"),"!")
  }
  if (NCOL(params)==1) {
    coef(object) <- params
    params <- as.matrix(params)
  }
  paramnames <- rownames(params)
  if (is.null(paramnames))
    pStop_(sQuote("params")," must have rownames")

  # create array to store weights per particle per block_list
  weights <- array(data = numeric(0), dim=c(nblocks,Np[1L]))
  loglik <- rep(NA,ntimes)

  for (nt in seq_len(ntimes)) {

    ## parameter autoregression
    foreach::foreach(par=sharedParNames,.combine=rbind) %do% {
      sharedPar <- params[paste0(par,1:U),,drop=FALSE]
      unit_mean <- apply(sharedPar,1,mean)
      overall_mean <- mean(unit_mean)
      reg <- sharedPar + spat_regression*(overall_mean - unit_mean)
      rownames(reg) <- paste0(par,1:U)
      reg
    } -> params_regressed
    if(!is.null(sharedParNames)){
      params[sharedParNamesExpanded,] <- params_regressed
    }
    ## parameter perturbation
    pmag <- cooling.fn(nt,nbpf)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation_spatPomp',params,pmag)
    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
    # note: params is on the estimation scale; tparams on the natural scale

     ## get initial states
      if (nt == 1L) {
        Xfilt <- rinit(object,params=tparams)
      }

    Xpred <- tryCatch(
      rprocess(
        object,
        x0=Xfilt,
        t0=times[nt],
        times=times[nt+1],
        params=tparams,
        .gnsi=gnsi
      ),
      error = function (e) {
        pStop_("process simulation error: ", conditionMessage(e))
      }
    )


    # this chunk of code is the same as for bpfilter
    max_log_d <- vector(mode = "numeric", length = nblocks)
    # For each  block, get each particle's weight
    for(i in seq(nblocks)){
      block <- block_list[[i]]
      log_vd <- tryCatch(
        vec_dmeasure(
          object,
          y=object@data[,nt,drop=FALSE],
          x=Xpred,
          units=block,
          times=times[nt+1],
          params=tparams,
          log=TRUE,
          .gnsi=gnsi
        ),
        error = function (e) {
          pStop_("error in calculation of weights: ",
               conditionMessage(e))
        }
      )
      log_d <- apply(log_vd[,,1,drop=FALSE], 2, function(x) sum(x))
      max_log_d[i] <- max(log_d)
      log_d <- log_d - max_log_d[i]
      weights[i,] <- exp(log_d)


    }
    gnsi <- FALSE

    ## resample for each block
    ## in future, this loop could be put into C if it is slow 
    for(i in seq_len(nblocks)){

      block = block_list[[i]]
      statenames_block = paste0(rep(object@unit_statenames,length(block)),
        rep(block,each=length(object@unit_statenames)))
      X_block = Xpred[statenames_block,,,drop = FALSE]

      paramnames_block = paste0(rep(estParNames,length(block)),
        rep(block,each=length(estParNames)))
      params_block <- params[paramnames_block,,drop=FALSE]

      xx <- tryCatch( #resampling with cross pollination
        .Call(
          "bpfilter_computations",
          x=X_block,
          params=params_block,
          Np=Np[nt+1],
          trackancestry=FALSE,
          doparRS=TRUE,
          weights=weights[i,]
        ),
        error = function (e) {
          pStop_(conditionMessage(e)) 
        }
      )

      Xfilt[statenames_block,] <- xx$states
      params[paramnames_block,] <- xx$params
    }
    
    
    log_weights = max_log_d + log(weights)
    loglik[nt] = sum(apply(log_weights,1,logmeanexp))

    if(nt==ntimes){
      mean_by_particles <- apply(params,1L,mean)
      coef(object,transform=TRUE) <- mean_by_particles
    }
  } ## end of main loop

  new(
    "ibpfd_spatPomp",
    object,
    block_list=block_list,
    Np=as.integer(Np),
    paramMatrix = params,
    cond.loglik=loglik,
    loglik=sum(loglik)
  )
}
