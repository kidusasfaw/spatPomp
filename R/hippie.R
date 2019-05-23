##' High-dimensional Iterated Perturbed Parameter Inference Engine (HIPPIE)
##'
##' An algorithm for estimating the parameters of a spatiotemporal partially-observed Markov process.
##' Running \code{hippie} causes the algorithm to perform a specified number of iterations of adapted simulations with parameter perturbation and parameter resamplings.
##' At each iteration, adapted simulations are performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' After cycling through the data, each island's weight is calculated and is used to rank the islands. The highest ranking islands are recycled into the next iteration.
##' This extra variability introduced through parameter perturbation effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##'
##' @name hippie
##' @rdname hippie
##' @include spatPomp_class.R
##' @family particle filter methods
##' @family \pkg{spatPomp} parameter estimation methods
##'
##' @inheritParams spatPomp
##' @inheritParams pfilter
##' @param Nhippie The number of iterations to perform.
##' @param islands The number of islands for the adapted simulations.
##' @param Np the number of particles inside each island.
##' @param rw.sd specification of the magnitude of the random-walk perturbations that will be applied to some or all model parameters.
##' Parameters that are to be estimated should have positive perturbations specified here.
##' The specification is given using the \code{\link{rw.sd}} function, which creates a list of unevaluated expressions.
##' The latter are evaluated in a context where the model time variable is defined (as \code{time}).
##' If \code{transform = TRUE}, then they are applied on the estimation scale.
##' @param cooling.type,cooling.fraction.50 specifications for the cooling schedule,
##' i.e., the manner and rate with which the intensity of the parameter perturbations is reduced with successive filtering iterations.
##' \code{cooling.type} specifies the nature of the cooling schedule.
##' See below (under \dQuote{Specifying the perturbations}) for more detail.
##'
##' @return
##' Upon successful completion, \code{hippie} returns an object of class
##' \sQuote{hippied.spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ extracts the point estimate }
##' }
##'
##' @section Specifying the perturbations:
##' The \code{rw.sd} function simply returns a list containing its arguments as unevaluated expressions.
##' These are then evaluated in a context containing the model \code{time} variable.  This allows for easy specification of the structure of the perturbations that are to be applied.
##' For example,
##' \preformatted{
##'     rw.sd(a=0.05, b=ifelse(0.2,time==time[1],0),
##'           c=ivp(0.2), d=ifelse(time==time[13],0.2,0),
##'           e=ivp(0.2,lag=13), f=ifelse(time<23,0.02,0))
##' }
##' results in perturbations of parameter \code{a} with s.d. 0.05 at every time step, while parameters \code{b} and \code{c} both get perturbations of s.d. 0.2 only before the first observation.
##' Parameters \code{d} and \code{e}, by contrast, get perturbations of s.d.  0.2 only before the thirteenth observation.
##' Finally, parameter \code{f} gets a random perturbation of size 0.02 before every observation falling before \eqn{t=23}.
##'
##' On the \eqn{m}-th HIPPIE iteration, prior to time-point \eqn{n}, the \eqn{d}-th parameter is given a random increment normally distributed with mean \eqn{0} and standard deviation \eqn{c_{m,n} \sigma_{d,n}}{c[m,n] sigma[d,n]}, where \eqn{c} is the cooling schedule and \eqn{\sigma}{sigma} is specified using \code{rw.sd}, as described above.
##' Let \eqn{N} be the length of the time series and \eqn{\alpha=}{alpha=}\code{cooling.fraction.50}.
##' Then, when \code{cooling.type="geometric"}, we have \deqn{c_{m,n}=\alpha^{\frac{n-1+(m-1)N}{50N}}.}{c[m,n]=alpha^((n-1+(m-1)N)/(50N)).}
##' When \code{cooling.type="hyperbolic"}, we have \deqn{c_{m,n}=\frac{s+1}{s+n+(m-1)N},}{c[m,n]=(s+1)/(s+n+(m-1)N),} where \eqn{s} satisfies \deqn{\frac{s+1}{s+50N}=\alpha.}{(s+1)/(s+50N)=alpha.}
##' Thus, in either case, the perturbations at the end of 50 HIPPIE iterations are a fraction \eqn{\alpha}{alpha} smaller than they are at first.
##'
NULL
rw.sd <- pomp2:::safecall


## define the hippied.spatPomp class
setClass(
  'hippied.spatPomp',
  contains='spatPomp',
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
  slots=c(
    log.island.weight="numeric",
    island.param="numeric",
    indices="vector",
    Np="integer",
    tol="numeric"
  ),
  prototype=prototype(
    log.island.weight=numeric(0),
    island.param=numeric(0),
    indices=integer(0),
    Np=as.integer(NA),
    tol=as.double(NA)
  )
)

hippie_pfilter.internal <- function (object, params, Np,
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
    log.island.weight=log.island.weight,
    island.param=params[,1],
    Np=as.integer(Np),
    tol=tol,
    indices=.indices
  )
}




hippie.internal <- function (object, islands, prop, Nhippie, start, Np, rw.sd, transform = FALSE,
                           cooling.type, cooling.fraction.50,
                           tol = (1e-18)^17, max.fail = Inf, island_bound = 10000,
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
  island.param.list = list()
  for(i in seq_len(islands)){
    .indices[[i]] <- integer(0)
    if (is.null(.paramMatrix)) {
      if (.ndone > 0) {               # call is from 'continue'
        island.param.list[[i]] <- coef(object)
        start <- coef(object)
      } else {                      # initial call
        island.param.list[[i]] <- start
      }
    } else {
      island.param.list[[i]] <- .paramMatrix[,i]
      start <- apply(.paramMatrix,1L,mean)
    }

    conv.rec <- array(dim=c(Nhippie+1,length(start)+2),
                      dimnames=list(iteration=seq.int(.ndone,.ndone+Nhippie),
                                    variable=c('loglik','nfail',names(start))))
    conv.rec[1L,] <- c(NA,NA,start)

    if (transform)
      island.param.list[[i]] <- partrans(object,island.param.list[[i]],dir="toEst",.gnsi=gnsi)
  }
  object <- as(object,"spatPomp")

  # iterate the filtering
  # cores <- parallel:::detectCores()-1
  doParallel::registerDoParallel(cores=NULL)
  mcopts <- list(set.seed=TRUE)
  mult.island.output <- list()

  for (n in seq_len(Nhippie)) {
    # begin multi-threaded
    print(n)
    weights <- vector(length=islands)
    param.swarm = array(dim=c(length(start),islands), dimnames=list(var = names(start), island = 1:islands))
    mult.island.output <- foreach::foreach(i=1:islands, .options.multicore=mcopts) %dopar%  {
      hippie_pfilter.internal(
        object=object,
        params=island.param.list[[i]],
        Np=Np,
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
    for(i in 1:islands){
      param.swarm[,i] <- mult.island.output[[i]]@island.param
      island.param.list[[i]] <- mult.island.output[[i]]@island.param
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

  pompUnload(object,verbose=verbose)

  new(
    "hippied.spatPomp",
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
  signature=signature(object="spatPomp"),
  definition = function (object, Nhippie = 1, islands, prop, start, Np,
                         rw.sd, transform = FALSE,
                         cooling.type = c("hyperbolic", "geometric"),
                         cooling.fraction.50,
                         tol = (1e-18)^17, max.fail = Inf, island_bound = 10000,
                         verbose = getOption("verbose"),...) {

    ep <- paste0("in ",sQuote("hippie"),": ")

    Nhippie <- as.integer(Nhippie)
    if (missing(start)) start <- coef(object)
    if (length(start)==0)
      stop(ep,sQuote("start")," must be specified if ",
           sQuote("coef(object)")," is NULL",call.=FALSE)
    if (is.null(names(start)))
      stop(ep,sQuote("start")," must be a named vector",
           call.=FALSE)
    ntimes <- length(time(object))
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
    if (missing(rw.sd))
      stop(ep,sQuote("rw.sd")," must be specified!",call.=FALSE)
    rw.sd <- pomp2:::perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
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
      rw.sd=rw.sd,
      transform=transform,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      tol=tol,
      max.fail=max.fail,
      island_bound = island_bound,
      verbose=verbose,
      ...
    )
  }
)
