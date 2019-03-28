##
## @references
## J. Park & E. L. Ionides (2018).
## A guided intermediate resampling particle filter for inference on high dimensional systems.
## Arxiv 1708.08543v2
##' @include spatpomp_class.R
##

setClass(
  "girfd_spatpomp",
  contains="spatpomp",
  slots=c(
    Ninter="numeric",
    Nguide="numeric",
    lookahead="numeric",
    h = "function",
    theta.to.v = "function",
    v.to.theta = "function",
    paramMatrix="array",
    eff.sample.size="array",
    cond.loglik="array",
    saved.states="array",
    Np="integer",
    tol="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    Ninter=as.double(NA),
    Nguide=as.double(NA),
    lookahead=as.double(NA),
    h = function(){},
    theta.to.v = function(){},
    v.to.theta = function(){},
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    eff.sample.size=array(data=numeric(0),dim=c(0,0)),
    cond.loglik=array(data=numeric(0),dim=c(0,0)),
    saved.states=array(data=numeric(0), dim = c(0,0)),
    Np=as.integer(NA),
    tol=as.double(NA),
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
    theta.to.v,
    v.to.theta,
    tol = 1e-17, max.fail = Inf,
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
        theta.to.v,
        v.to.theta,
        tol = tol, max.fail = Inf,
        save.states = FALSE,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp2:::pStop("girf",conditionMessage(e))
    )
  }
)

##' @name girf-girfd_spatpomp
##' @aliases girf,girfd_spatpomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(data="girfd_spatpomp"),
  function (data,
            Np,
            tol,
            Ninter,
            Nguide,
            lookahead,
            h,
            theta.to.v,
            v.to.theta,
            ...,
            verbose = getOption("verbose", FALSE)
            ) {
    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol
    if (missing(Ninter)) Ninter <- data@Ninter
    if (missing(Nguide)) Nguide <- data@Nguide
    if (missing(lookahead)) lookahead <- data@lookahead
    if (missing(h)) h <- data@h
    if (missing(theta.to.v)) theta.to.v <- data@theta.to.v
    if (missing(v.to.theta)) v.to.theta <- data@v.to.theta

    girf(as(data,"spatpomp"),
         Np=Np,
         tol=tol,
         Ninter=Ninter,
         Nguide=Nguide,
         lookahead=lookahead,
         h=h,
         theta.to.v=theta.to.v,
         v.to.theta=v.to.theta,
         ...,
         verbose=verbose)

  }
)

girf.internal <- function (object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        h,
        theta.to.v,
        v.to.theta,
        tol = 1e-17, max.fail = Inf,
        save.states = FALSE,
        cooling, cooling.m,
        .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)
  ep <- paste0("in ",sQuote("girf"),": ")

  if (pomp2:::undefined(object@rprocess) || pomp2:::undefined(object@dmeasure))
    pomp2:::pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  save.states <- as.logical(save.states)

  params <- coef(object)
  params.matrix <- matrix(params,nrow=length(params), ncol = Np[1])
  rownames(params.matrix) <- names(params)
  times <- time(object,t0=TRUE)
  t0 <- times[1]
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
  if (save.states) {
    xparticles <- setNames(vector(mode="list",length=ntimes),time(object))
  }

  cond.loglik <- array(0, dim = c(ntimes, Ninter))
  eff.sample.size <- array(0, dim = c(ntimes, Ninter))

  # initialize filter guide function
  filter_guide_fun <- array(1, dim = Np[1])
  for (nt in 0:(ntimes-1)) { ## main loop
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    # four-dimensional array: nvars by nguide by ntimes by nreps
    Xg = array(0, dim=c(length(statenames), Nguide, lookahead_steps, Np[1]), dimnames = list(nvars = statenames, ng = NULL, lookahead = 1:lookahead_steps, nreps = NULL))
    ## for each particle get K guide particles, and fill in sample variance over K for each (lookahead value - unit - particle) combination
    fcst_samp_var <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
    for (p in 1:Np[1]){
      # find this particle's initialization and repeat in Nguide times
      xp = matrix(x[,p], nrow = nrow(x), ncol = Nguide, dimnames = list(nvars = statenames, ng = NULL))
      # get all the guides for this particles
      Xg[,,,p] <- rprocess(object, x0=xp, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
               params=params,.gnsi=gnsi)
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for(l in 1:lookahead_steps){
          hXg = apply(X=Xg[snames,,l,p, drop = FALSE], MARGIN = c(2,3,4), FUN = h, obj=object)
          fcst_samp_var[u, l, p] = var(hXg)
        }
      }
    }
    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      # get prediction simulations
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
                    params=params,.gnsi=gnsi)

      # X is now a nvars by nreps by 1 array
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp2::flow(object, x0=X.start, t0=tt[s+1], params=params.matrix, times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)])
      } else {
        skel <- X
      }

      # create measurement variance at skeleton matrix
      meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for (l in 1:lookahead_steps){
          hskel <- apply(X=skel[snames,,l, drop = FALSE], MARGIN = c(2,3), FUN = h, obj = object)
          meas_var_skel[u,l,] <- apply(X=hskel, MARGIN = c(1,2), FUN = theta.to.v, obj = object)
        }
      }

      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        for(l in 1:lookahead_steps){
          fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
                                      FUN = function(x) x*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1]))
        }
      }

      mom_match_param <- array(0, dim = c(length(params), length(object@units), lookahead_steps, Np[1]), dimnames = list(params = names(params), lookahead = NULL, J = NULL))
      inflated_var <- meas_var_skel + fcst_var_upd
      mom_match_param = apply(X=inflated_var, MARGIN=c(1,2,3), FUN = v.to.theta, obj = object)
      # guide functions as product (so base case is 1)
      guide_fun = vector(mode = "numeric", length = Np[1]) + 1

      for(l in 1:lookahead_steps){
        dmeas_weights <- tryCatch(
          vec_dmeasure(
            object,
            y=object@data[,nt+l,drop=FALSE],
            x=skel[,,l,drop = FALSE],
            times=times[nt+1+l],
            params=mom_match_param[,,l,],
            log=FALSE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        resamp_weights <- apply(dmeas_weights[,,1,drop=FALSE], 2, function(x) prod(x))
        guide_fun = guide_fun*resamp_weights
      }

      guide_fun[guide_fun < tol^(lookahead*length(object@units))] <- tol^(lookahead*length(object@units))
      s_not_1_weights <- guide_fun/filter_guide_fun
      if (!(s==1 & nt!=0)){
        weights <- s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        weights <- tryCatch(
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=params,
            log=FALSE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        weights <- as.numeric(weights)*s_not_1_weights
      }

      xx <- tryCatch(
        .Call('girf_computations',
              x=X,
              params=params,
              Np=Np[nt+1],
              predmean=FALSE,
              predvar=FALSE,
              filtmean=FALSE,
              trackancestry=FALSE,
              doparRS=FALSE,
              weights=weights,
              gps=guide_fun,
              fsv=fcst_samp_var,
              tol=tol
              ),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      all.fail <- xx$fail
      eff.sample.size[nt+1, s] <- xx$ess
      cond.loglik[nt+1, s] <- xx$loglik
      x <- xx$states
      filter_guide_fun <- xx$filterguides
      params <- xx$params[,1]
      fcst_samp_var <- xx$newfsv
    }
  }
  new(
    "girfd_spatpomp",
    object,
    Ninter=Ninter,
    Nguide=Nguide,
    lookahead=lookahead,
    h = h,
    theta.to.v = theta.to.v,
    v.to.theta = v.to.theta,
    paramMatrix=params.matrix,
    eff.sample.size=eff.sample.size,
    cond.loglik=cond.loglik,
    saved.states=x,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}
