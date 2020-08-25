setClass(
  "bpfilterd_spatPomp",
  contains="spatPomp",
  slots=c(
    partition="list",
    Np="integer",
    cond.loglik="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    partition = list(),
    Np=as.integer(NA),
    cond.loglik=as.double(NA),
    loglik=as.double(NA)
  )
)

##' Block particle filter (BPF)
##'
##' An algorithm used to estimate the filter distribution of a spatiotemporal partially-observed Markov process (spatPomp)
##' Running \code{bpfilter} causes the algorithm to split the spatial units into different partitions so that each spatial
##' unit belongs to one partition. After the particles are propagated, resampling of the particles occurs
##' within each partition independently based on sampled weights within the partition. Each partition samples only the spatial
##' components within the partition which allows for cross-pollination of particles where the highest weighted
##' components of each particle are more likely to be resampled and get combined with resampled components of other particles.
##' By using local particle filters and resampling with a smaller subset of dimensions, it tries to avert the curse of dimensionality so that
##' the resampling does not result in particle depletion with one particle representing the complex filter distribution.
##'
##' @name bpfilter-spatPomp
##' @aliases bpfilter,spatPomp-method
##' @rdname bpfilter
##' @include spatPomp_class.R generics.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @param object A \code{spatPomp} object
##' @param params A parameter set for the spatiotemporal POMP. If missing, \code{bpfilter} will attempt to run using \code{coef(object)}
##' @param Np The number of particles used for the simulations. If missing, \code{bpfilter} will attempt to run using \code{ncol(params)}
##' @param block_size The number of spatial units per block.
##' @param partition List that can specifies a partition of the spatial units. Each partition element called a \code{block} and is
##' an integer vector of neighboring units.
##'
##' @examples
##' # Create a simulation of a BM using default parameter set
##' b <- bm(U=6, N=10)
##'
##' # Run BPF with the specified number of particles and number of units per block
##' bpfilterd.b1 <- bpfilter(b, Np = 100, block_size = 2)
##'
##' # Run BPF with the specified number of particles and partition. This specification
##' is exactly equivalent to the previous example
##' bpfilterd.b2 <- bpfilter(b, Np = 20, partition = list(c(1,2), c(3,4), c(5,6)))
##'
##' # Get a likelihood estimate
##' logLik(bpfilterd.b2)
##'
##' @return
##' Upon successful completion, \code{bpfilter} returns an object of class
##' \sQuote{bpfilterd_spatPomp}.
##'
##' @section Details:
##' Only one of \code{block_size} or \code{partition} should be specified.
##' If both or neither is provided, an error is triggered.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields an estimate of the log-likelihood of the data under the model. }
##' }
##'
##' @export
setMethod(
  "bpfilter",
  signature=signature(object="spatPomp"),
  function (object, Np, block_size, partition, params) {
    ep = paste0("in ",sQuote("bpfilter"),": ")

    if (missing(params)) params <- coef(object)

    if(missing(partition) && missing(block_size))
      stop(ep,sQuote("partition"), " or ", sQuote("block_size"), " must be specified to the call",call.=FALSE)

    if (!missing(partition) & !missing(block_size)){
      stop(ep,"Exactly one of ",sQuote("block_size"), " and ", sQuote("partition"), " should be provided, but not both.",call.=FALSE)
    }

    if (missing(Np)) {
      if (is.matrix(params)) {
        Np <- ncol(params)
      } else {
        stop(ep,sQuote("Np")," must be specified",call.=FALSE)
      }
    }

    if (missing(partition)){
      if(block_size > length(unit_names(object))){
        stop(ep,sQuote("block_size"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      all_units = seq_len(length(unit_names(object)))
      nblocks = round(length(all_units)/block_size)
      partition = split(all_units, sort(all_units %% nblocks))
    }
    partition <- lapply(partition, as.integer)

    bpfilter.internal(
     object=object,
     Np=Np,
     partition=partition,
     params=params)
  }
)
bpfilter.internal <- function (object, Np, partition, params, .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("bpfilter"),": ")
  object <- as(object,"spatPomp")
  pompLoad(object)
  gnsi <- as.logical(.gnsi)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))
  nblocks <- length(partition)

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
    one.par <- TRUE
    coef(object) <- params
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

  # create array to store weights per particle per partition
  weights <- array(data = numeric(0), dim=c(nblocks,Np[1L]))
  loglik <- rep(NA,ntimes)

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

    # For each  partition, get each particle's weight
    for(i in seq(nblocks)){
      block <- partition[[i]]
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

    ## resample for each partition
    for(i in seq_len(nblocks)){
      block = partition[[i]]
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
          weights=weights[i,]
        ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      x[statenames,] <- xx$states
      params <- xx$params
    }
    log_weights = max_log_d + log(weights)
    loglik[nt] = sum(apply(log_weights,1,logmeanexp))
  } ## end of main loop
  new(
    "bpfilterd_spatPomp",
    object,
    partition=partition,
    Np=as.integer(Np),
    cond.loglik=loglik,
    loglik=sum(loglik)
  )
}
