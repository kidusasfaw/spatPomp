## HIPPIE algorithm functions

rw.sd <- pomp2:::safecall


## define the hippied.pomp class
setClass(
  'hippied.spatpomp',
  contains='spatpomp',
  slots=c(
    Nhippie = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    transform = 'logical',
    conv.rec = 'matrix',
    resamp.frac = 'numeric',
    paramMatrix = 'array'
  )
)

setClass(
  "hippie.as.pfilterd.pomp",
  contains="spatpomp",
  slots=c(
    #loc.comb.pred.weights="numeric",
    #cond.densities="array",
    log.island.weight="numeric",
    #paramMatrix="array",
    island.param="numeric",
    indices="vector",
    #eff.sample.size="numeric",
    #cond.loglik="numeric",
    Np="integer",
    tol="numeric"
    #nfail="integer",
    #loglik="numeric"
  ),
  prototype=prototype(
    #loc.comb.pred.weights = numeric(0),
    #cond.densities=array(data=numeric(0),dim=c(0,0,0)),
    log.island.weight=numeric(0),
    #paramMatrix=array(data=numeric(0),dim=c(0,0)),
    island.param=numeric(0),
    indices=integer(0),
    #eff.sample.size=numeric(0),
    #cond.loglik=numeric(0),
    saved.states=list(),
    saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA)
    #nfail=as.integer(NA),
    #loglik=as.double(NA)
  )
)

hippie_pfilter.internal <- function (object, params, Np, nbhd,
                          hippieiter, cooling.fn, rw.sd,
                          tol = (1e-18)^17, max.fail = Inf,
                          transform, .indices = integer(0), verbose,
                          .getnativesymbolinfo = TRUE) {
  ep <- paste0("in ",sQuote("hippie_pfilter.internal"),": ")
  gnsi <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)
  verbose <- as.logical(verbose)
  hippieiter <- as.integer(hippieiter)
  Np <- as.integer(Np)

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    stop(ep,sQuote(".indices"),
         " has improper length",call.=FALSE)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  units <- unit(object)
  nunits <- length(units)

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  params <- array(data = params, dim = c(length(params),Np[1L]), dimnames = list(variable=names(params),rep=NULL))

  # create array to store weights across time
  cond.densities <- array(data = numeric(0), dim=c(nunits,Np[1L],ntimes))
  dimnames(cond.densities) <- list(unit = 1:nunits, rep = 1:Np[1L], time = 1:ntimes)
  log.island.weight <- 0
  for (nt in seq_len(ntimes)) {

    ## perturb parameters
    pmag <- cooling.fn(nt,hippieiter)$alpha*rw.sd[,nt]
    params <- .Call('randwalk_perturbation',params,pmag,PACKAGE = 'pomp2')

    if (transform)
      tparams <- pomp2::partrans(object,params,dir="fromEst",
                          .gnsi=gnsi)

    if (nt == 1L) {
      ## get initial states
      # x <- init.state(object,nsim=Np[1L],params=if (transform) tparams else params)
      x <- rinit(object,params=tparams)
    }

    ## advance the state variables according to the process model
    X <- tryCatch(
      rprocess(
        object,
        xstart=x,
        times=times[c(nt,nt+1)],
        params=if (transform) tparams else params,
        offset=1,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,"process simulation error: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    ## determine the weights. returns weights which is a nunits by Np by ntimes array
    weights <- tryCatch(
      vec_dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=X,
        times=times[nt+1],
        params=if (transform) tparams else params,
        log=FALSE,
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,"error in calculation of weights: ",
             conditionMessage(e),call.=FALSE)
      }
    )

    weights[is.na(weights)] <- tol
    resamp_weights <- apply(weights[,,1,drop=FALSE], 2, function(x) prod(x))
    # if any particle's resampling weight is zero divide out it's weight vector by the smallest component
    if(all(resamp_weights == 0)) resamp_weights <- rep(tol, Np[1L])
    log.island.weight <- log.island.weight + log((1/Np[1L])*sum(resamp_weights))
    gnsi <- FALSE

    xx <- tryCatch(
      .Call(
        hippie_computations,
        x=X,
        params=params,
        Np=Np[nt+1],
        rw_sd=numeric(0),
        predmean=FALSE,
        predvar=FALSE,
        filtmean=FALSE,
        trackancestry=do_ta,
        onepar=FALSE,
        weights=resamp_weights
        #tol=tol
      ),
      error = function (e) {
        stop(ep,"hippie pfilter computations error: ",conditionMessage(e),call.=FALSE)
      }
    )

    if (do_ta) {
      .indices <- .indices[xx$ancestry]
    }

    x <- xx$states
    params <- xx$params

    if (nt == ntimes) {
      if (any(resamp_weights>0)) {
        coef(object,transform=transform) <- apply(params,1L,weighted.mean,w=resamp_weights)
      } else {
        warning(ep,"filtering failure at last filter iteration, using unweighted mean for ",
                sQuote("coef"),call.=FALSE)
        coef(object,transform=transform) <- apply(params,1L,mean)
      }
    }

    if (verbose && (nt%%5==0))
      cat("hippie pfilter timestep",nt,"of",ntimes,"finished\n")
  }


  new(
    "hippie.as.pfilterd.pomp",
    object,
    log.island.weight=log.island.weight,
    island.param=params[,1],
    Np=as.integer(Np),
    tol=tol,
    indices=.indices
  )
}




hippie.internal <- function (object, islands, prop, Nhippie, start, Np, nbhd, rw.sd, transform = FALSE,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17, max.fail = Inf,
                           verbose = FALSE, .ndone = 0L,
                           .indices = vector(mode="list", length = islands),
                           .paramMatrix = NULL,
                           .getnativesymbolinfo = TRUE, ...) {

  ep <- paste0("in ",sQuote("hippie"),": ")

  pompLoad(object,verbose=verbose)

  transform <- as.logical(transform)
  verbose <- as.logical(verbose)
  gnsi <- as.logical(.getnativesymbolinfo)
  Np <- c(Np,Np[1L])


  if (Nhippie <= 0)
    stop(ep,sQuote("Nhippie")," must be a positive integer",call.=FALSE)

  if (islands <= 0)
    stop(ep,sQuote("islands")," must be a positive integer",call.=FALSE)

  cooling.fn <- pomp2:::mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )
  # initialize the parameter for each island
  # paramMatrixList = list()
  island.param.list = list()
  for(i in seq_len(islands)){
    .indices[[i]] <- integer(0)
    if (is.null(.paramMatrix)) {
      if (.ndone > 0) {               # call is from 'continue'
        # paramMatrixList[[i]] <- object@paramMatrix
        island.param.list[[i]] <- coef(object)
        #start <- apply(paramMatrixList[[i]],1L,mean)
        start <- coef(object)
      } else {                      # initial call
        # paramMatrixList[[i]] <- array(data=start,dim=c(length(start),Np[1L]), dimnames=list(variable=names(start),rep=NULL))
        island.param.list[[i]] <- start
      }
    } else {
      # paramMatrixList[[i]] <- .paramMatrix
      island.param.list[[i]] <- .paramMatrix[,i]
      # start <- apply(paramMatrixList[[i]],1L,mean)
      start <- apply(.paramMatrix,1L,mean)
    }

    conv.rec <- array(dim=c(Nhippie+1,length(start)+2),
                      dimnames=list(iteration=seq.int(.ndone,.ndone+Nhippie),
                                    variable=c('loglik','nfail',names(start))))
    conv.rec[1L,] <- c(NA,NA,start)

    if (transform)
      # paramMatrixList[[i]] <- partrans(object,paramMatrixList[[i]],dir="toEstimationScale", .getnativesymbolinfo=gnsi)
      island.param.list[[i]] <- partrans(object,island.param.list[[i]],dir="toEst",.gnsi=gnsi)
  }
  object <- as(object,"spatpomp")

  # iterate the filtering
  require(doParallel)
  require(parallel)
  # cores <- parallel:::detectCores()-1
  registerDoParallel(cores=NULL)
  mcopts <- list(set.seed=TRUE)
  mult.island.output <- list()

  for (n in seq_len(Nhippie)) {
    # param.swarm = array(dim=c(length(start),islands), dimnames=list(var = names(start), island = 1:islands))
    # begin multi-threaded
    mult.island.output <- foreach(i=1:islands, .options.multicore=mcopts) %dopar%  {
      hippie_pfilter.internal(
        object=object,
        # params = paramMatrixList[[i]]
        params=island.param.list[[i]],
        Np=Np,
        nbhd=nbhd,
        hippieiter=.ndone+n,
        cooling.fn=cooling.fn,
        rw.sd=rw.sd,
        tol=tol,
        max.fail=max.fail,
        transform=transform,
        .indices=.indices[[i]],
        verbose=verbose
      )
    }
    # end multi-threaded
    # begin single-threaded
    # for(i in 1:islands){
    #   mult.island.output[[i]] <- hippie_pfilter.internal(
    #     object=object,
    #     #params=paramMatrixList[[i]],
    #     params=island.param.list[[i]]
    #     Np=Np,
    #     nbhd=nbhd,
    #     hippieiter=.ndone+n,
    #     cooling.fn=cooling.fn,
    #     rw.sd=rw.sd,
    #     tol=tol,
    #     max.fail=max.fail,
    #     transform=transform,
    #     .indices=.indices[[i]],
    #     verbose=verbose
    #   )
    # }
    # end single threaded
    gnsi <- FALSE
    fails <- 0
    weights <- vector(length=islands)
    param.swarm = array(dim=c(length(start),islands), dimnames=list(var = names(start), island = 1:islands))
    for(i in 1:islands){
      param.swarm[,i] <- mult.island.output[[i]]@island.param
      island.param.list[[i]] <- mult.island.output[[i]]@island.param
      # fails <- fails + mult.island.output[[i]]@nfail
      .indices[[i]] <- mult.island.output[[i]]@indices
      weights[i] <- mult.island.output[[i]]@log.island.weight
    }

    coef(object, transform = transform) <- apply(param.swarm,1,mean)

    conv.rec[n+1,-c(1,2)] <- coef(object)
    conv.rec[n,c(1,2)] <- c(0,0)

    # top p quantile weights
    top_indices <- which(weights > quantile(weights, p = 1-prop))
    assignments <- rep(0, islands)
    assignments <- suppressWarnings(top_indices + assignments) # uses recycling to ensure equitable distribution
    assignments <- sample(assignments)


    for(i in 1:islands) {
      # in cases where all weights are equal and extremely low, this will not happen
      if(length(assignments) > 0) island.param.list[[i]] <- island.param.list[[assignments[i]]]
    }

    if (n != Nhippie){
      rm(param.swarm)
      rm(assignments)
      rm(top_indices)
      rm(weights)
    }

    if (verbose) cat("hippie iteration",n,"of",Nhippie,"completed\n")
  }

  # parameter swarm to be outputted
  if (transform)
    param.swarm <- partrans(object,param.swarm,dir="fromEst",
                                .gnsi=gnsi)
  #coef(object, transform=transform) <- apply(final.param.swarm,1L,mean) # weighted mean hard because weights unstable.

  pompUnload(object,verbose=verbose)

  new(
    "hippied.spatpomp",
    object,
    Nhippie=Nhippie,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    transform=transform,
    conv.rec=conv.rec,
    paramMatrix=param.swarm
  )
}

setMethod(
  "hippie",
  signature=signature(object="spatpomp"),
  definition = function (object, Nhippie = 1, islands, prop, start, Np,
                         nbhd, rw.sd, transform = FALSE,
                         cooling.type = c("hyperbolic", "geometric"),
                         cooling.fraction.50,
                         tol = (1e-18)^17, max.fail = Inf,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("hippie"),": ")

    Nhippie <- as.integer(Nhippie)
    print("before coef object")
    if (missing(start)) start <- coef(object)
    if (length(start)==0)
      stop(ep,sQuote("start")," must be specified if ",
           sQuote("coef(object)")," is NULL",call.=FALSE)
    if (is.null(names(start)))
      stop(ep,sQuote("start")," must be a named vector",
           call.=FALSE)
    print("before time object")
    ntimes <- length(time(object))
    print("after time object")
    if (missing(Np))
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    else if (is.function(Np)) {
      Np <- tryCatch(
        vapply(seq_len(ntimes),Np,numeric(1)),
        error = function (e) {
          stop(ep,"if ",sQuote("Np"),
               " is a function, it must return a single positive integer",
               call.=FALSE)
        }
      )
    } else if (!is.numeric(Np))
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
    if(missing(islands))
      stop(ep,"number of islands, ",
           sQuote("islands"),", must be specified!",call.=FALSE)
    if(missing(prop))
      stop(ep,"top proportion of islands (by weight) to be resampled, ",
           sQuote("prop"),", must be specified!",call.=FALSE)
    if(missing(nbhd)){
      nbhd <- function(object, time, unit) {
        nunits = length(unit(object))
        ntimes = length(time(object)) #length(time(object))
        nbhd_matrix = array(0, dim = c(nunits, ntimes))
        # B_{d,n} = entire past
        if(time > 1) {
          nbhd_matrix[1:unit, 1:(time-1)] = 1
        }
        return(nbhd_matrix)
      }
    }
    print("before perturbation kernel")
    if (missing(rw.sd))
      stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
    # rw.sd <- pomp2:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
    rw.sd <- pomp2:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
    print("after perturbation kernel")
    cooling.type <- match.arg(cooling.type)

    cooling.fraction.50 <- as.numeric(cooling.fraction.50)
    if (cooling.fraction.50 <= 0 || cooling.fraction.50 > 1)
      stop(ep,sQuote("cooling.fraction.50"),
           " must be in (0,1]",call.=FALSE)

    hippie.internal(
      object=object,
      islands=islands,
      prop=prop,
      Nhippie=Nhippie,
      start=start,
      Np=Np,
      nbhd=nbhd,
      rw.sd=rw.sd,
      transform=transform,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      ...
    )
  }
)

setMethod(
  "hippie",
  signature=signature(object="hippied.spatpomp"),
  definition = function (object, Nhippie=1, start, Np,
                         rw.sd, transform, cooling.type, cooling.fraction.50,
                         tol, ...) {

    if (missing(Nhippie)) Nhippie <- object@Nhippie
    if (missing(start)) start <- coef(object)
    if (missing(rw.sd)) rw.sd <- object@rw.sd
    if (missing(transform)) transform <- object@transform
    if (missing(cooling.type)) cooling.type <- object@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- object@cooling.fraction.50

    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol

    f <- selectMethod("hippie","pomp")

    f(object,Nhippie=Nhippie,start=start,Np=Np,rw.sd=rw.sd,transform=transform,
      cooling.type=cooling.type,cooling.fraction.50=cooling.fraction.50,
      tol=tol,...)
  }
)

setMethod(
  'continue',
  signature=signature(object='hippied.spatpomp'),
  definition = function (object, Nhippie = 1, ...) {

    ndone <- object@Nhippie

    f <- selectMethod("hippie","hippied.pomp")
    obj <- f(object=object,Nhippie=Nhippie,.ndone=ndone,...)

    object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1L,c('loglik','nfail')]
    obj@conv.rec <- rbind(
      object@conv.rec,
      obj@conv.rec[-1L,colnames(object@conv.rec)]
    )
    names(dimnames(obj@conv.rec)) <- c("iteration","variable")
    obj@Nhippie <- as.integer(ndone+Nhippie)

    obj
  }
)
