## independent island filtering codes
setClass(
  "island.spatpomp",
  contains="pomp",
  slots=c(
    loc.comb.pred.weights="array",
    cond.densities="array",
    pred.mean="array",
    pred.var="array",
    filter.mean="array",
    filter.traj="array",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    saved.states="list",
    saved.params="list",
    Np="integer",
    tol="numeric",
    nfail="integer",
    loglik="numeric"
  ),
  prototype=prototype(
    loc.comb.pred.weights=array(data=numeric(0),dim=c(0,0)),
    cond.densities=array(data=numeric(0),dim=c(0,0,0)),
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    pred.var=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    filter.traj=array(data=numeric(0),dim=c(0,0,0)),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    saved.states=list(),
    saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=as.integer(NA),
    loglik=as.double(NA)
  )
)
setClass(
  "iifd.spatpomp",
  contains="pomp",
  slots=c(
    loc.comb.filter.weights="array",
    pred.mean="list",
    pred.var="list",
    filter.mean="list",
    filter.traj="list",
    paramMatrix="array",
    indices="vector",
    eff.sample.size="list",
    cond.loglik="array",
    saved.states="list",
    saved.params="list",
    Np="integer",
    tol="numeric",
    nfail="list",
    loglik="numeric"
  ),
  prototype=prototype(
    loc.comb.filter.weights=array(data=numeric(0),dim=c(0,0,0,0)),
    pred.mean=list(),
    pred.var=list(),
    filter.mean=list(),
    filter.traj=list(),
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    eff.sample.size=list(),
    cond.loglik=array(data=numeric(0),dim=c(0,0)),
    saved.states=list(),
    saved.params=list(),
    Np=as.integer(NA),
    tol=as.double(NA),
    nfail=list(),
    loglik=as.double(NA)
  )
)
iif.internal <- function (object, params, Np,
  nbhd,
  tol, max.fail,
  pred.mean = FALSE,
  pred.var = FALSE,
  filter.mean = FALSE,
  filter.traj = FALSE,
  cooling, cooling.m,
  verbose = FALSE,
  save.states = FALSE,
  save.params = FALSE,
  .getnativesymbolinfo = TRUE) {

  ep <- paste0("in ",sQuote("iif"),": ")
  object <- as(object,"spatpomp")
  pompLoad(object,verbose=verbose)
  gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
  pred.mean <- as.logical(pred.mean)
  pred.var <- as.logical(pred.var)
  filter.mean <- as.logical(filter.mean)
  filter.traj <- as.logical(filter.traj)
  verbose <- as.logical(verbose)
  save.states <- as.logical(save.states)
  save.params <- as.logical(save.params)

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  if (missing(tol))
    stop(ep,sQuote("tol")," must be specified",call.=FALSE)

  one.par <- FALSE
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  units <- unit(object)
  nunits <- length(units)

  if (missing(Np)) {
    if (is.matrix(params)) {
      Np <- ncol(params)
    } else {
      stop(ep,sQuote("Np")," must be specified",call.=FALSE)
    }
  }
  if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        stop(ep,"if ",sQuote("Np")," is a function, ",
          "it must return a single positive integer",call.=FALSE)
      }
    )
  }
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

  if (NCOL(params)==1) {        # there is only one parameter vector
    one.par <- TRUE
    coef(object) <- params     # set params slot to the parameters
    params <- as.matrix(params)
  }

  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  ## init.state returns an nvars by nsim matrix
  init.x <- init.state(object,params=params,nsim=Np[1L])
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x

  ## set up storage for saving samples from filtering distributions
  if (save.states | filter.traj) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (save.params) {
    pparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  } else {
    pparticles <- list()
  }
  if (filter.traj) {
    pedigree <- vector(mode="list",length=ntimes+1)
  }

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  ## set up storage for prediction means, variances, etc.
  if (pred.mean) {
    pred.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (pred.var) {
    pred.v <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    pred.v <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.mean) {
    filt.m <- matrix(
      data=0,
      nrow=nvars,
      ncol=ntimes,
      dimnames=list(
        variable=statenames,
        time=time(object))
    )
  } else {
    filt.m <- array(data=numeric(0),dim=c(0,0))
  }

  if (filter.traj) {
    filt.t <- array(
      data=0,
      dim=c(nvars,1,ntimes+1),
      dimnames=list(
        variable=statenames,
        rep=1,
        time=times)
    )
  } else {
    filt.t <- array(data=numeric(0),dim=c(0,0,0))
  }

  # create array to store weights across time
  cond.densities <- array(data = numeric(0), dim=c(nunits,Np[1L],ntimes))
  dimnames(cond.densities) <- list(unit = 1:nunits, rep = 1:Np[1L], time = 1:ntimes)
  #preweighted <- vector("list", length = ntimes)
  #resampling_weights <- vector("list", length = ntimes)
  #postweighted <- vector("list", length = ntimes)

  # before rprocess, let us convert covariates to wide-format
}
