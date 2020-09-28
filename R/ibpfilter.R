##' Iterated block particle filter (IBPF)
##'
##' An implementation of a parameter estimation algorithm combining
##' GIRF with IF2, proposed by Park and Ionides (2019), following the pseudocode in Asfaw, Ionides and King (2019).
##'
##' @name ibpfilter
##' @rdname ibpfilter
##' @include spatPomp_class.R generics.R spatPomp.R bpfilter.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @inheritParams spatPomp
##' @inheritParams pomp::mif2
##'
##' @param Nbpf the number of iterations of perturbed BPF.
##'
##' @return
##' Upon successful completion, \code{ibpfilter} returns an object of class
##' \sQuote{ibpfilterd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ gives the Monte Carlo estimate of the maximum likelihood. }
##' }
##'

NULL

rw.sd <- pomp:::safecall

setClass(
  "ibpfilterd_spatPomp",
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
  "ibpfilter",
  function (data, ...)
    standardGeneric("ibpfilter")
)

setMethod(
  "ibpfilter",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("ibpfilter","data")
  }
)

setMethod(
  "ibpfilter",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("ibpfilter",data)
  }
)

##' @name ibpfilter-spatPomp
##' @aliases ibpfilter,spatPomp-method
##' @rdname ibpfilter
##' @export
setMethod(
  "ibpfilter",
  signature=signature(data="spatPomp"),
  definition=function (data,Nbpf,Np,rw.sd,cooling.type,cooling.fraction.50,
                       block_size, block_list,spat_regression,
                       tol = 1e-300, max.fail = Inf,save.states = FALSE,
                       ..., verbose = getOption("verbose", FALSE)) {
    if(missing(block_list) && missing(block_size))
      stop(ep,sQuote("block_list"), " or ", sQuote("block_size"), " must be specified to the call",call.=FALSE)

    if (!missing(block_list) & !missing(block_size)){
      stop(ep,"Exactly one of ",sQuote("block_size"), " and ", sQuote("block_list"), " should be provided, but not both.",call.=FALSE)
    }

    if (missing(spat_regression))
      stop(ep, sQuote("spat_regression"), " should be provided",call.=FALSE)

    if (missing(block_list)){
      if(block_size > length(unit_names(data))){
        stop(ep,sQuote("block_size"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      all_units = seq_len(length(unit_names(data)))
      nblocks = round(length(all_units)/block_size)
      block_list = split(all_units, sort(all_units %% nblocks))
    }
    block_list <- lapply(block_list, as.integer)


    if (missing(Np)) {
      if (is.matrix(params)) {
        Np <- ncol(params)
      } else {
        stop(ep,sQuote("Np")," must be specified",call.=FALSE)
      }
    }

    tryCatch(
      ibpfilter_internal(data,Nbpf,spat_regression,Np,rw.sd,cooling.type,cooling.fraction.50,
                     block_list,tol = tol,max.fail = Inf, save.states = FALSE,...,verbose=verbose),
      error = function (e) pomp:::pStop("ibpfilter",conditionMessage(e))
    )
  }
)

ibpfilter_internal <- function (object,Nbpf,spat_regression,Np,rw.sd,cooling.type,cooling.fraction.50,
                            block_list,tol, max.fail = Inf,save.states = FALSE,...,
                            .ndone = 0L, .indices = integer(0),.paramMatrix = NULL,.gnsi = TRUE, verbose = FALSE) {

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
  if (pomp:::undefined(object@rprocess) || pomp:::undefined(object@dunit_measure))
    pStop_(paste(sQuote(c("rprocess","dunit_measure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Nbpf) != 1 || !is.numeric(Nbpf) || !is.finite(Nbpf) || Nbpf < 1)
    pStop_(sQuote("Nbpf")," must be a positive integer.")
  Nbpf <- as.integer(Nbpf)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))
  nunits <- length(unit_names(object))

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
      pWarn("igirf","Np[k] ignored for k > ",sQuote("length(time(object))"),".")
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
  rw.sd <- pomp:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- pomp:::mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]*nunits),
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

  ## iterate the filtering
  for (n in seq_len(Nbpf)) {
    b <- ibpfilter_bpfilter(object=object,block_list=block_list,params=paramMatrix,
                    spat_regression=spat_regression,Np=Np,bpfiter=.ndone+n,cooling.fn=cooling.fn,rw.sd=rw.sd,tol=tol,max.fail=max.fail,
                    verbose=verbose,.indices=.indices, .gnsi=gnsi)
    gnsi <- FALSE

    paramMatrix <- b@paramMatrix
    traces[n+1,-1] <- coef(b)
    traces[n+1,1] <- b@loglik
    .indices <- .indices

    if (verbose) cat("ibpfilter iteration",n,"of",Nbpf,"completed\n")

  }

  b@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
                            .gnsi=gnsi)

  new(
    "ibpfilterd_spatPomp",
    b,
    Nbpf=Nbpf,
    spat_regression=spat_regression,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

ibpfilter_bpfilter <- function (object,block_list,params,spat_regression,Np,bpfiter,cooling.fn,rw.sd, tol,max.fail = Inf,
                                verbose,.indices = integer(0), .gnsi = TRUE) {
  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  bpfiter <- as.integer(bpfiter)
  Np <- as.integer(Np)
  ep <- paste0("in ",sQuote("ibpfilter"),": ")
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small positive number.")

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))
  nblocks <- length(block_list)
  loglik <- rep(NA,ntimes)

  if (length(Np)==1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np)!=(ntimes+1))
    stop(ep,sQuote("Np")," must have length 1 or length ",ntimes+1,call.=FALSE)
  if (any(Np<=0))
    stop(ep,"number of particles, ",sQuote("Np"),", must always be positive",call.=FALSE)
  if (!is.numeric(Np))
    stop(ep,sQuote("Np")," must be a number, a vector of numbers, or a function",call.=FALSE)
  Np <- as.integer(Np)
  print(Np)
  print(dim(params))
  if (is.matrix(params)) {
    if (!all(Np*nunits==ncol(params)))
      stop(ep,"when ",sQuote("params")," is provided as a matrix, do not specify ",
           sQuote("Np"),"!",call.=FALSE)
  }
  if (NCOL(params)==1) {
    one.par <- TRUE
    coef(object) <- params
    params <- as.matrix(params)
  }
  paramnames <- rownames(params)
  npars <- length(paramnames)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  # create array to store weights per particle per block_list
  weights <- array(data = numeric(0), dim=c(nblocks,Np[1L]))
  loglik <- rep(NA,ntimes)

  # 1  2  3  4  5  6  7  8  9 10 11 12 to 1  4  7 10  2  5  8 11  3  6  9 12
  # used to go from npars by (U * J) where the first U columns belong to j=1
  # to go to npars by (J * U) where the first J columns belong to u=1
  u_by_np <- matrix(seq.int(nunits*Np[1L]),nrow=nunits,ncol=Np[1L])
  grouped_by_u_ix <- as.numeric(t(u_by_np))
  grouped_by_np_ix <- as.numeric(u_by_np)
  # to get a dim(P) by U matrix of parameters for each unit,
  # average for each unit the J particles that represent that unit's parameters
  params_by_unit <- t(rowsum(t(params), as.integer(gl(ncol(params), Np[1L], ncol(params))))) / Np[1L]
  spat_weight_mat <- diag(spat_regression, ncol = nunits, nrow = nunits)
  spat_weight_mat[upper.tri(spat_weight_mat) | lower.tri(spat_weight_mat)] <- 1-spat_regression
  # perform the spatial regression to get initial parameters for this bpfilter iteration
  params_new <- params_by_unit %*% spat_weight_mat
  rownames(params_new) <- paramnames
  for (nt in seq_len(ntimes)) {
    if (nt == 1L) {
      ## perturb parameters
      pmag <- cooling.fn(nt,bpfiter)$alpha*rw.sd[,nt]
      params <- .Call('randwalk_perturbation',params_new,pmag,PACKAGE = 'pomp')
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
      tparamnames <- rownames(tparams)
      # the first particle is now the first U columns
      tparams <- matrix(rep(tparams, Np[1L]) , nrow = nrow(tparams), byrow=FALSE)
      dimnames(tparams) <- list(tparamnames, NULL)
      # get tparams by particle by averaging over the U parameters for each particle
      tparams_by_particle <- t(rowsum(t(tparams), as.integer(gl(ncol(tparams), nunits, ncol(tparams))))) / nunits
      dimnames(tparams_by_particle) <- list(tparamnames, NULL)
      # get initial states
      X <- rinit(object,params=tparams_by_particle)
      xnames <- rownames(X)
    }
    else{
      pmag <- cooling.fn(nt,bpfiter)$alpha*rw.sd[,nt]
      params <- partrans(object,tparams,dir="toEst",.gnsi=gnsi)
      params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp')
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
    }

    # get tparams by particle by averaging over the U parameters for each particle
    tparams_by_particle <- t(rowsum(t(tparams), as.integer(gl(ncol(tparams), nunits, ncol(tparams))))) / nunits
    dimnames(tparams_by_particle) <- list(tparamnames, NULL)

    X <- tryCatch(
      rprocess(
        object,
        x0=X,
        t0=times[nt],
        times=times[nt+1],
        params=tparams_by_particle,
        .gnsi=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    # prepare tparams_by_particle to be passed to vec_dmeasure
    tparams_by_units_by_particle <- tparams_by_particle
    dim(tparams_by_units_by_particle) <- c(npars,nunits,Np[1L])
    dimnames(tparams_by_units_by_particle) <- list(tparamnames, NULL, NULL)

    max_log_d <- vector(mode = "numeric", length = nblocks)
    # For each  block, get each particle's weight
    for(i in seq(nblocks)){
      block <- block_list[[i]]
      log_vd <- tryCatch(
        vec_dmeasure(
          object,
          y=object@data[,nt,drop=FALSE],
          x=X,
          units=block,
          times=times[nt+1],
          params=tparams_by_units_by_particle[,block,,drop=FALSE],
          log=TRUE,
          .gnsi=gnsi
        ),
        error = function (e) {
          stop(ep,"error in calculation of weights: ",
               conditionMessage(e),call.=FALSE)
        }
      )
      log_d <- apply(log_vd[,,1,drop=FALSE], 2, function(x) sum(x))
      max_log_d[i] <- max(log_d)
      log_d <- log_d - max_log_d[i]
      weights[i,] <- exp(log_d)
    }
    gnsi <- FALSE

    ## resample for each block
    for(i in seq_len(nblocks)){
      block = block_list[[i]]
      stacked_tparams_by_particle <- tparams_by_units_by_particle
      dim(stacked_tparams_by_particle) <- c(npars*nunits, Np[1L])
      dimnames(stacked_tparams_by_particle) <- list(rep(paramnames,nunits),NULL)
      temp_tparams <- stacked_tparams_by_particle[(1 + (i-1)*npars*length(block)):(i*npars*length(block)),drop=FALSE]
      rownames(temp_tparams) <- tparamnames
      us = object@unit_statenames
      statenames = paste0(rep(us,length(block)),rep(block,each=length(us)))
      tempX = X[statenames,,,drop = FALSE]

      xx <- tryCatch( #resampling with cross pollination
        .Call(
          "bpfilter_computations",
          x=tempX,
          params=temp_tparams,
          Np=Np[nt+1],
          trackancestry=FALSE,
          doparRS=TRUE,
          weights=weights[i,]
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      X[statenames,] <- xx$states
      stacked_tparams_by_particle[(1 + (i-1)*npars*length(block)):(i*npars*length(block)),] <- xx$params
    }
    log_weights = max_log_d + log(weights)
    loglik[nt] = sum(apply(log_weights,1,logmeanexp))
    tparams <- stacked_tparams_by_particle
    dim(tparams) <- c(npars,nunits*Np[1L])
    dimnames(tparams) <- list(tparamnames,NULL)

    if(nt==ntimes){
      coef(object,transform=FALSE) <- apply(tparams,1L,mean)
      params <- partrans(object,tparams,dir="toEst",.gnsi=gnsi)
    }
  } ## end of main loop
  new(
    "bpfilterd_spatPomp",
    object,
    block_list=block_list,
    Np=as.integer(Np),
    # go back to the first J particles belonging to u=1 case
    paramMatrix = params[,grouped_by_u_ix,drop=FALSE],
    cond.loglik=loglik,
    loglik=sum(loglik)
  )
}
