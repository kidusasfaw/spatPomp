##
## @references
## J. Park & E. L. Ionides (2018).
## A guided intermediate resampling particle filter for inference on high dimensional systems.
## Arxiv 1708.08543v2
##

setClass(
  "girfd_pomp",
  contains="pomp",
  slots=c(
    S="integer",
    K="integer",
  # remaining slots are from pfilter
    pred.mean="array",
    pred.var="array",
    filter.mean="array",
    filter.traj="array",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    saved.states="list",
    Np="integer",
    tol="numeric",
    nfail="integer",
    loglik="numeric"
  ),
  prototype=prototype(
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    pred.var=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    saved.states=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=as.integer(NA),
    loglik=as.double(NA)
  )
)

setGeneric(
  "girf",
  function (data, ...)
    standardGeneric("girf")
)

setMethod(
  "girf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("girf","data")
  }
)

setMethod(
  "girf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("girf",data)
  }
)

##' @name girf-spatpomp
##' @aliases girf,spatpomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(data="spatpomp"),
  definition=function (
    data,
    Np,
    Ninter,
    lookahead,
    Nguide,
    h,
    theta_to_v,
    v_to_theta,
    tol = 1e-17, max.fail = Inf,
    pred.mean = FALSE,
    pred.var = FALSE,
    filter.mean = FALSE,
    filter.traj = FALSE,
    save.states = FALSE,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      girf.internal(
        data,
        Np,
        Ninter,
        lookahead,
        Nguide,
        h,
        theta_to_v,
        v_to_theta,
        tol = 1e-17, max.fail = Inf,
        pred.mean = FALSE,
        pred.var = FALSE,
        filter.mean = FALSE,
        filter.traj = FALSE,
        save.states = FALSE,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp2:::pStop("girf",conditionMessage(e))
    )
  }
)

##' @name girf-girfd_pomp
##' @aliases girf,girfd_pomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(data="girfd_pomp"),
  function (data, Np, tol,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    girf(as(data,"pomp"),Np=Np,tol=tol,
      ...,verbose=verbose)

  }
)

girf.internal <- function (object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        h,
        theta_to_v,
        v_to_theta,
        tol = 1e-17, max.fail = Inf,
        pred.mean = FALSE,
        pred.var = FALSE,
        filter.mean = FALSE,
        filter.traj = FALSE,
        save.states = FALSE,
        cooling, cooling.m,
        .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)

  if (pomp2:::undefined(object@rprocess) || pomp2:::undefined(object@dmeasure))
    pomp2:::pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  save.states <- as.logical(save.states)

  params <- coef(object)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np) || is.null(Np)) {
    pomp2:::pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pomp2:::pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pomp2:::pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pomp2:::pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pomp2:::pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pomp2:::pStop_(sQuote("tol")," should be a small positive number.")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states || filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- array(data=numeric(1),dim=c(nvars,ntimes),
      dimnames=list(variable=statenames,time=time(object)))
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(data=numeric(1),dim=c(nvars,1,ntimes+1),
      dimnames=list(variable=statenames,rep=1,time=times))
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  for (nt in 1:ntimes) { ## main loop
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt],to=times[nt+1],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    # four-dimensional array: nvars by nguide by ntimes by nreps
    Xg = array(0, dim=c(length(statenames), Nguide, lookahead_steps, Np[1]), dimnames = list(nvars = statenames, ng = NULL, lookahead = 1:lookahead_steps, nreps = NULL))
    ## for each particle get K guide particles, and fill in sample variance over K for each (lookahead value - unit - particle) combination
    fcst_samp_var <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
    for (p in 1:Np[1]){
      # find this particle's initialization and repeat in Nguide times
      xp = matrix(x[,p], nrow = nrow(x), ncol = Nguide, dimnames = dimnames(x))
      # get all the guides for this particles
      Xg[,,,p] <- rprocess(object, xstart=xp, times=times[nt:(nt+lookahead)],
               params=params,offset=1L,.gnsi=gnsi)
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for(l in 1:lookahead_steps){
          hXg = apply(X=Xg[snames,,l,p, drop = FALSE], MARGIN = c(2,3,4), FUN = h)
          fcst_samp_var[u, l, p] = var(hXg)
        }
      }
    }
    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      # get prediction simulations
      X <- rprocess(object,xstart=x,times=c(tt[s], tt[s+1]),
                    params=params,offset=1L,.gnsi=gnsi)
      # X is now a nvars by nreps by 1 array
      X.start <- X[,,1]
      skel <- pomp2::flow(object, xstart=X.start, tstart = tt[s+1], params=params, times = c(tt[s+1], times[(nt+1):(nt + lookahead_steps)]))
      # create measurement variance at skeleton matrix
      meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for (l in 1:lookahead_steps){
          hskel <- apply(X=skel[snames,,l, drop = FALSE], MARGIN = c(2,3), FUN = h)
          meas_var_skel[u,l,] <- apply(X=hskel, MARGIN = c(1,2), FUN = theta_to_v)
        }
      }
      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        for(l in 1:lookahead_steps){
          fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
                                      FUN = function(x) x*(times[nt+l] - tt[s+1])/(times[nt+l] - times[nt]))
        }
      }
      mom_match_param <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      inflated_var <- meas_var_skel + fcst_var_upd
      theta <- v_to_theta(inflated_var)

      # X.skel <- X[,,]

      # for(l in 1:lookahead_steps){
      #   skel <- skeleton(object, x=X, times = c(tt[s], times[nt + lookahead_steps]), params=params)
      #   return(skel)
      # }
    }
  }
}
#     ## advance the state variables according to the process model
#     X <- rprocess(object,xstart=x,times=times[c(nt,nt+1)],params=params,
#       offset=1L,.gnsi=gnsi)
#
#     if (pred.var) { ## check for nonfinite state variables and parameters
#       problem.indices <- unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
#       nprob <- length(problem.indices)
#       if (nprob > 0)
#         pomp2:::pStop_("non-finite state variable",ngettext(nprob,"","s"),": ",
#           paste(rownames(X)[problem.indices],collapse=', '))
#     }
#
#     ## determine the weights
#     weights <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
#       times=times[nt+1],params=params,log=FALSE,.gnsi=gnsi)
#
#     if (!all(is.finite(weights))) {
#       first <- which(!is.finite(weights))[1L]
#       datvals <- object@data[,nt]
#       weight <- weights[first]
#       states <- X[,first,1L]
#       msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,
#         datvals,states,params)
#       pomp2:::pStop_(msg)
#     }
#
#     gnsi <- FALSE
#
#     ## compute prediction mean, prediction variance, filtering mean,
#     ## effective sample size, log-likelihood
#     ## also do resampling if filtering has not failed
#     xx <- .Call(P_girf_computations,x=X,params=params,Np=Np[nt+1],
#       predmean=pred.mean,predvar=pred.var,filtmean=filter.mean,
#       trackancestry=filter.traj,doparRS=FALSE,weights=weights,tol=tol)
#
#     all.fail <- xx$fail
#     loglik[nt] <- xx$loglik
#     eff.sample.size[nt] <- xx$ess
#
#     x <- xx$states
#     params <- xx$params[,1L]
#
#     if (pred.mean) pred.m[,nt] <- xx$pm
#     if (pred.var) pred.v[,nt] <- xx$pv
#     if (filter.mean) filt.m[,nt] <- xx$fm
#     if (filter.traj) pedigree[[nt]] <- xx$ancestry
#
#     if (all.fail) { ## all particles are lost
#       nfail <- nfail+1
#       if (verbose) message("filtering failure at time t = ",times[nt+1])
#       if (nfail>max.fail) pomp2:::pStop_("too many filtering failures")
#     }
#
#     if (save.states || filter.traj) {
#       xparticles[[nt]] <- x
#       dimnames(xparticles[[nt]]) <- setNames(dimnames(xparticles[[nt]]),
#         c("variable","rep"))
#     }
#
#     if (verbose && (nt%%5 == 0))
#       cat("girf timestep",nt,"of",ntimes,"finished\n")
#
#   } ## end of main loop
#
#   if (filter.traj) { ## select a single trajectory
#     b <- sample.int(n=length(weights),size=1L,replace=TRUE)
#     filt.t[,1L,ntimes+1] <- xparticles[[ntimes]][,b]
#     for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
#       b <- pedigree[[nt+1]][b]
#       filt.t[,1L,nt+1] <- xparticles[[nt]][,b]
#     }
#     if (times[2L] > times[1L]) {
#       b <- pedigree[[1L]][b]
#       filt.t[,1L,1L] <- init.x[,b]
#     } else {
#       filt.t <- filt.t[,,-1L,drop=FALSE]
#     }
#   }
#
#   if (!save.states) xparticles <- list()
#
#   if (nfail>0)
#     pWarn("girf",nfail," filtering failure",ngettext(nfail,"","s")," occurred.")
#
#   new(
#     "girfd_pomp",
#     object,
#     pred.mean=pred.m,
#     pred.var=pred.v,
#     filter.mean=filt.m,
#     filter.traj=filt.t,
#     paramMatrix=array(data=numeric(0),dim=c(0,0)),
#     eff.sample.size=eff.sample.size,
#     cond.loglik=loglik,
#     saved.states=xparticles,
#     Np=as.integer(Np),
#     tol=tol,
#     nfail=as.integer(nfail),
#     loglik=sum(loglik)
#   )
# }

nonfinite_dmeasure_error <- function (time, lik, datvals, states, params) {
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",preserve.width="common")
  paste0(
    sQuote("dmeasure")," returns non-finite value.\n",
    "likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}
