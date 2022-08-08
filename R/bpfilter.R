##' Block particle filter (BPF)
##'
##' An implementation of the block particle filter algorithm of Rebeschini and van Handel (2015), which is used to estimate the filter distribution
##' of a spatiotemporal partially-observed Markov process.
##' \code{bpfilter} requires a partition of the spatial units which can be provided by either the \code{block_size} or the \code{block_list} argument.
##' The elements of the partition are called blocks. We perform resampling for each block independently based on sample weights within the block.
##' Each resampled block only contains latent states for the spatial components within the block which allows for a ``cross-pollination" of
##' particles where the highest weighted segments of each particle are more likely to be resampled and get combined with resampled components of
##' other particles. The method mitigates the curse of dimensionality by resampling locally.
##'
##' @name bpfilter
##' @rdname bpfilter
##' @include spatPomp_class.R
##' @family particle filter methods
##' @inheritParams abf
##'
##' @param block_size The number of spatial units per block. If this is provided, the method subdivides units approximately evenly
##' into blocks with size \code{block_size}.
##' @param block_list List that specifies an exact partition of the spatial units. Each partition element, or block, is
##' an integer vector of neighboring units.
##' @param save_states logical. If True, the state-vector for each particle and
##' block is saved.
##' @param \dots If a \code{params} argument is specified, \code{bpfilter} will estimate the likelihood at that parameter set instead of at \code{coef(object)}.
##'
##' @examples
##' # Create a simulation of a Brownian motion
##' b <- bm(U=4, N=2)
##'
##' # Run BPF with the specified number of units per block
##' bpfilterd_b1 <- bpfilter(b, Np = 10, block_size = 2)
##'
##' # Run BPF with the specified partition
##' bpfilterd_b2 <- bpfilter(b,
##'                          Np = 10,
##'                          block_list = list(c(1,2),c(3,4)) )
##'
##' # Get a likelihood estimate
##' logLik(bpfilterd_b2)
##'
##' @return Upon successful completion, \code{bpfilter()} returns an object of class
##' \sQuote{bpfilterd_spatPomp} containing the algorithmic parameters used to run \code{bpfilter()}
##' and the estimated likelihood.
##'
##' @section Details:
##' Only one of \code{block_size} or \code{block_list} should be specified.
##' If both or neither is provided, an error is triggered.
##'
##' @references \rebeschini2015
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields an estimate of the log-likelihood of the data under the model. }
##' }
##'

setClass(
  "bpfilterd_spatPomp",
  contains="spatPomp",
  slots=c(
    block_list="list",
    Np="integer",
    paramMatrix="array",
    cond.loglik="numeric",
    block.cond.loglik="array",
    loglik="numeric",
    saved.states="list"
  ),
  prototype=prototype(
    block_list = list(),
    Np=as.integer(NA),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    cond.loglik=as.double(NA),
    block.cond.loglik=array(data=numeric(0),dim=c(0,0)),
    loglik=as.double(NA),
    saved.states=list()
  )
)

setGeneric("bpfilter",  function (object, ...) standardGeneric("bpfilter"))

##' @name bpfilter-spatPomp
##' @aliases bpfilter,spatPomp-method
##' @rdname bpfilter
##' @export
setMethod(
  "bpfilter",
  signature=signature(object="spatPomp"),
  function (object, Np, block_size, block_list, save_states, ..., verbose=getOption("verbose", FALSE)) {
    ep = paste0("in ",sQuote("bpfilter"),": ")

    if(missing(save_states)) save_states <- FALSE

    if(missing(block_list) && missing(block_size))
      stop(ep,sQuote("block_list"), " or ", sQuote("block_size"), " must be specified to the call",call.=FALSE)

    if (!missing(block_list) & !missing(block_size)){
      stop(ep,"Exactly one of ",sQuote("block_size"), " and ", sQuote("block_list"), " should be provided, but not both.",call.=FALSE)
    }

    if (missing(Np)) {
        stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    }

    if (missing(block_list)){
      if(block_size > length(unit_names(object))){
        stop(ep,sQuote("block_size"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      all_units = seq_len(length(unit_names(object)))
      nblocks = round(length(all_units)/block_size)
      block_list = split(all_units, sort(all_units %% nblocks))
    }
    block_list <- lapply(block_list, as.integer)

    bpfilter.internal(
     object=object,
     Np=Np,
     block_list=block_list,
     save_states=save_states,
     ...,
     verbose=verbose)
  }
)

bpfilter.internal <- function (object, Np, block_list, save_states, ..., verbose, .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("bpfilter"),": ")
  verbose <- as.logical(verbose)
  p_object <- pomp(object,...)
  save_states <- as.logical(save_states)
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
                unit_obsnames = object@unit_obsnames,
                unit_accumvars = object@unit_accumvars)
  params <- coef(object)
  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))
  gnsi <- as.logical(.gnsi)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nblocks <- length(block_list)

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
  if (NCOL(params)==1) {
    coef(object) <- params
    params <- as.matrix(params)
  }
  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  ## returns an nvars by nsim matrix
  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  x <- init.x

  # create array to store weights per particle per block_list
  weights <- array(data = numeric(0), dim=c(nblocks,Np[1L]))
  loglik <- rep(NA,ntimes)
  block.loglik <- matrix(NA,nblocks,ntimes)
  if (save_states) {
    saved.states <- setNames(vector(mode="list",length=ntimes), time(object))
  } else {
    saved.states <- list()
  }


  for (nt in seq_len(ntimes)) { ## main loop
    ## advance the state variables according to the process model
    max_log_d <- vector(mode = "numeric", length = nblocks)
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
          params=params,
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
      us = object@unit_statenames
      statenames = paste0(rep(us,length(block)),rep(block,each=length(us)))
      tempX = X[statenames,,,drop = FALSE]
      xx <- tryCatch( #resampling with cross pollination
        .Call(
          "bpfilter_computations",
          x=tempX,
          params=params,
          Np=Np[nt+1],
          trackancestry=FALSE,
          doparRS=FALSE,
          weights=weights[i,]
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      x[statenames,] <- xx$states
      params <- xx$params
    }
    if (save_states) saved.states[[nt]] <- x
    log_weights = max_log_d + log(weights)
    block_log_weights <- apply(log_weights,1,logmeanexp)
    loglik[nt] = sum(block_log_weights)
    block.loglik[,nt] <- block_log_weights
  } ## end of main loop
  new(
    "bpfilterd_spatPomp",
    object,
    block_list=block_list,
    Np=as.integer(Np),
    cond.loglik=loglik,
    block.cond.loglik=block.loglik,
    loglik=sum(loglik),
    saved.states=saved.states
  )
}
