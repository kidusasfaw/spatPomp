setClass(
  "bpfilterd_spatPomp",
  contains="spatPomp",
  slots=c(
    partitions_list="list",
    Np="integer",
    cond.loglik="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    partitions_list = list(),
    Np=as.integer(NA),
    cond.loglik=as.double(NA),
    loglik=as.double(NA)
  )
)

##' @name bpfilter-spatPomp
##' @aliases bpfilter,spatPomp-method
##' @rdname bpfilter
##' @export
setMethod(
  "bpfilter",
  signature=signature(object="spatPomp"),
  function (object, Np, num_partitions, partitions_list, params) {
    ep = paste0("in ",sQuote("bpfilter"),": ")
    if (missing(params)) params <- coef(object)
    if(missing(partitions_list) && missing(num_partitions))
      stop(ep,sQuote("partitions_list"), " or ", sQuote("num_partitions"), " must be specified to the call",call.=FALSE)
    if (missing(Np)) {
      if (is.matrix(params)) {
        Np <- ncol(params)
      } else {
        stop(ep,sQuote("Np")," must be specified",call.=FALSE)
      }
    }
    if (missing(partitions_list)){
      if(num_partitions > length(spat_units(object))){
        stop(ep,sQuote("num_partitions"), " cannot be greater than the number of spatial units",call.=FALSE)
      }
      partition_size = floor(length(spat_units(object))/num_partitions)
      partitions_list = split(spat_units(object), ceiling(seq_along(spat_units(object))/partition_size))
    }
    bpfilter.internal(
     object=object,
     Np=Np,
     partitions_list=partitions_list,
     params=params)
  }
)
bpfilter.internal <- function (object, Np, partitions_list, params, .gnsi = TRUE) {
  ep <- paste0("in ",sQuote("bpfilter"),": ")
  object <- as(object,"spatPomp")
  pompLoad(object)
  gnsi <- as.logical(.gnsi)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(spat_units(object))

  num_partitions <- length(partitions_list)

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
  weights <- array(data = numeric(0), dim=c(num_partitions,Np[1L]))
  loglik <- rep(NA,ntimes)

  for (nt in seq_len(ntimes)) { ## main loop
    ## advance the state variables according to the process model
    #cat("time \n", nt, "\n")
    max_log_d <- vector(mode = "numeric", length = num_partitions)
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

    #cat(" \n X \n")
    #print(X)


    # For each  partition, get each particle's weight
    for(i in seq(num_partitions)){
      partition <- partitions_list[[i]]
      # print(partition)
      log_vd <- tryCatch(
        vec_dmeasure(
          object,
          y=object@data[,nt,drop=FALSE],
          x=X,
          units=match(partition, spat_units(object)),
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
      # print("done with for loop")
      # print(dim(log_vd))
      # print(log_vd)
      log_d <- apply(log_vd[,,1,drop=FALSE], 2, function(x) sum(x))
      max_log_d[i] <- max(log_d)
      log_d <- log_d - max_log_d[i]
      weights[i,] <- exp(log_d)
    }
    gnsi <- FALSE
    #cat("\nweights \n")
    #print(weights)

    ## resample for each partition
    for(i in seq(num_partitions)){
      partition = partitions_list[[i]]
      statenames = paste0(paste0(rep(object@unit_statenames, each = length(partition))), match(partition, spat_units(object)))
      # cat("statenames \n", statenames, "\n")
      # cat("dim(X) \n", dim(X))
      tempX = X[statenames,,,drop = FALSE]
      # paste0("tempX \n", tempX, "\n")
      # print(Np)
      # paste0("Np \n", Np[nt+1], "\n")
      # print(weights)
      # paste0("weights \n", weights, "\n")

      xx <- tryCatch(
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
      # cat("statenames \n", statenames, "\n")
      # cat("dim(xx$states) \n", dim(xx$states), "\n")
      # cat("xx$states \n", xx$states, "\n")


      x[statenames,] <- xx$states

      params <- xx$params
    }
    #cat("\nnew x \n")
    #print(x)
    log_weights = max_log_d + log(weights)
    loglik[nt] = logmeanexp(apply(log_weights,2,sum))
  } ## end of main loop
  new(
    "bpfilterd_spatPomp",
    object,
    partitions_list=partitions_list,
    Np=as.integer(Np),
    cond.loglik=loglik,
    loglik=sum(loglik)
  )
}
