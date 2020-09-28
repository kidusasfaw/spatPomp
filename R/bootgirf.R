##' Guided intermediate resampling filter (GIRF) using a bootstrap guide function.
##'
##' An implementation of the algorithm of Park and Ionides (2020),
##' this function is under development, and later will be combined with girf(). In bootgirf.R (this file), the pseudo-simulations are obtained by adding difference in simulation residuals at two target times to the skeleton simulations.
##'
##' @name bootgirf
##' @rdname bootgirf
##' @include spatPomp_class.R generics.R spatPomp.R
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##'
##'
##' @inheritParams spatPomp
##' @param object A \code{spatPomp} object.
##' @param params A parameter set for the spatiotemporal POMP.
##' @param Np The number of Monte Carlo particles to be used.
##' @param Ninter the number of intermediate resampling timepoints.
##' @param lookahead The number of future observations included in the guide function.
##' @param Nguide The number of simulations used to estimate state process uncertainty for each particle.
##' @param tol If the guide functions become too small (beyond floating-point precision limits), we set them to this value.
##'
##' @examples
##' # Create a simulation of a BM using default parameter set
##' b <- bm(U=3, N=10)
##'
##' # Run bootstrap-GIRF
##' girfd.b <- bootgirf(b,
##'                 Np = 100,
##'                 Ninter = length(unit_names(b)),
##'                 lookahead = 1,
##'                 Nguide = 50
##' )
##' # Get the likelihood estimate from GIRF
##' logLik(girfd.b)
##'
##' # Compare with the likelihood estimate from particle filter
##' pfd.b <- pfilter(b, Np = 500)
##' logLik(pfd.b)
##' @return
##' Upon successful completion, \code{bootgirf} returns an object of class
##' \sQuote{girfd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields an unbiased estimate of the log-likelihood of the data under the model. }
##' }
##'
##' @references
##' \park2020
##'
##' \asfaw2019
NULL


setClass(
  "girfd_spatPomp",
  contains="spatPomp",
  slots=c(
    Ninter="numeric",
    Nguide="numeric",
    lookahead="numeric",
    cond.loglik="array",
    Np="integer",
    tol="numeric",
    loglik="numeric",
    paramMatrix="array"
  ),
  prototype=prototype(
    Ninter=as.double(NA),
    Nguide=as.double(NA),
    lookahead=as.double(NA),
    cond.loglik=array(data=numeric(0),dim=c(0,0)),
    Np=as.integer(NA),
    tol=as.double(NA),
    loglik=as.double(NA),
    paramMatrix=array(data=numeric(0),dim=c(0,0))
  )
)

setGeneric(
  "bootgirf",
  function (object, ...)
    standardGeneric("bootgirf")
)

setMethod(
  "bootgirf",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("bootgirf","data")
  }
)

setMethod(
  "bootgirf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("bootgirf",object)
  }
)

##' @name bootgirf-spatPomp
##' @aliases bootgirf,spatPomp-method
##' @rdname bootgirf
##' @export
setMethod(
  "bootgirf",
  signature=signature(object="spatPomp"),
  definition=function (
    object,
    Np,
    Ninter,
    lookahead,
    Nguide,
    params,
    tol,
    ...) {

    if (missing(params)) params <- coef(object)
    if (missing(tol)) tol <- 1e-300
    if (missing(Ninter)) Ninter <- length(unit_names(object))

    tryCatch(
      bootgirf.internal(
        object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        params,
        tol,
        ...
      ),
      error = function (e) pomp:::pStop("bootgirf",conditionMessage(e))
    )
  }
)

##' @name bootgirf-girfd_spatPomp
##' @aliases bootgirf,girfd_spatPomp-method
##' @rdname bootgirf
##' @export
setMethod(
  "bootgirf",
  signature=signature(object="girfd_spatPomp"),
  function (object,
            Np,
            Ninter,
            lookahead,
            Nguide,
            params,
            tol,
            ...
            ) {
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    if (missing(Ninter)) Ninter <- object@Ninter
    if (missing(Nguide)) Nguide <- object@Nguide
    if (missing(lookahead)) lookahead <- object@lookahead
    if (missing(params)) params <- coef(object)

    bootgirf(as(object,"spatPomp"),
         Np=Np,
         Ninter=Ninter,
         lookahead=lookahead,
         Nguide=Nguide,
         params = params,
         tol=tol,
         ...)

  }
)

bootgirf.internal <- function (object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        params,
        tol,
        ...,
        .gnsi = TRUE) {

  verbose <- FALSE
  ep <- paste0("in ",sQuote("bootgirf"),": ")

  if (pomp:::undefined(object@rprocess) || pomp:::undefined(object@dmeasure))
    pomp:::pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)

  params.matrix <- matrix(params,nrow=length(params), ncol = Np[1])
  rownames(params.matrix) <- names(params)
  times <- time(object,t0=TRUE)
  t0 <- times[1]
  ntimes <- length(times)-1
  U <- length(unit_names(object))
  if (missing(Np) || is.null(Np)) {
    pomp:::pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pomp:::pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pomp:::pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pomp:::pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pomp:::pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pomp:::pStop_(sQuote("tol")," should be a small positive number.")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  statenames <- rownames(init.x)
  nvars <- nrow(init.x)
  x <- init.x
  znames <- object@accumvars
  cond.loglik <- array(0, dim = c(ntimes, Ninter))
  # initialize filter guide function
  log_filter_guide_fun <- array(0, dim = Np[1]) 
  for (nt in 0:(ntimes-1)) { ## main loop
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    # get a matrix with nguides times nreps columns to propagate using rprocess
    x_with_guides <- x[,rep(1:Np[1], each=Nguide)]
    guidesim_times <- c(sapply(1:lookahead_steps, function(bb) seq(from=times[nt+bb],to=times[nt+bb+1],length.out=Ninter+1)[-1])) # times at which guide simulations will be recorded
    guidesim_index <- 1:Np[1] # the index for guide simulations (to be updated each time resampling occurs)
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=guidesim_times, params=params,.gnsi=gnsi)
    Xskel <- tryCatch( # skeleton 
          pomp::flow(object,
                     x0=x,
                     t0=times[nt+1],
                     params=params.matrix,
                     times = guidesim_times,
                     ...),
          error = function (e) {
            pomp::flow(object,
                       x0=x,
                       t0=times[nt+1],
                       params=params.matrix,
                       times = guidesim_times,
                       method = 'adams')
          }
        )
    resids <- Xg - Xskel[,rep(1:Np[1], each=Nguide),] # residuals

    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      # get prediction simulations
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
                    params=params,.gnsi=gnsi)
      # X is now a nvars by nreps by 1 array
      if(s>1 && length(znames)>0){
        x.znames <- x[znames,]; dim(x.znames) <- c(dim(x.znames),1)
        X[znames,,] <- X[znames,,,drop=FALSE] + x.znames
      }

      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- tryCatch(
          pomp::flow(object,
                     x0=X.start,
                     t0=tt[s+1],
                     params=params.matrix,
                     times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],
                     ...),
          error = function (e) {
            #stop(ep,conditionMessage(e),call.=FALSE) # nocov
            pomp::flow(object,
                       x0=X.start,
                       t0=tt[s+1],
                       params=params.matrix,
                       times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],
                       method = 'adams')
          }
        )
        if(length(znames) > 0){
          skel.start <- skel[,,1]
          X.start.znames <- X.start[znames,]
          skel.start.znames <- skel.start[znames,]
          skel.end.znames <- X.start.znames + skel.start.znames
          skel[znames,,1] <- skel.end.znames
        }
      } else {
        skel <- X
      }

      # guide functions as product (so base case is 1)
      log_guide_fun = vector(mode = "numeric", length = Np[1])

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead_steps <= 0) discount_denom_init = object@t0
        else discount_denom_init = times[nt+1+l - lookahead_steps]
        discount_factor = 1 - (times[nt+1+l] - tt[s+1])/max(times[nt+1+l] - discount_denom_init, 2*(times[nt+2]-times[nt+1])) ## the denominator is at least twice the observation interval, to ensure that the discount factor does not become too small for L=1 and small s (which can lead to very uninformative guide function.

        # construct pseudo-simulations by adding simulated noise terms (residuals) to the skeletons
        pseudosims <- skel[,rep(1:Np[1], each=Nguide),l,drop=FALSE] + resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),Ninter*l,drop=FALSE] - resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),s,drop=FALSE]
        log_dmeas_weights <- tryCatch(
          (vec_dmeasure(
            object,
            y=object@data[,nt+l,drop=FALSE],
            x=pseudosims,
            times=times[nt+1+l],
            params=params,
            log=TRUE, 
            .gnsi=gnsi
          )),
          error = function (e) {
            stop(ep,"error in calculation of log_dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        ldw <- array(log_dmeas_weights, c(U,Nguide,Np[1])) # log_dmeas_weights is an array with dim U*(Np*Nguide)*1. Reorder it as U*Nguide*Np
        log_fcst_lik <- colSums(log(apply(exp(ldw),c(1,3),sum)/Nguide)) # average dmeas (natural scale) over Nguide sims, then take log, and then sum over 1:U (for each particle)
        log_resamp_weights <- log_fcst_lik*discount_factor
        log_guide_fun = log_guide_fun + log_resamp_weights
      }
      log_s_not_1_weights <- log_guide_fun - log_filter_guide_fun
      if (!(s==1 & nt!=0)){
        log_weights <- log_s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        log_meas_weights <- tryCatch(
          (dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=params,
            log=TRUE, 
            .gnsi=gnsi
          )),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        log_weights <- as.numeric(log_meas_weights) + log_s_not_1_weights
      }
      max_log_weights <- max(log_weights, na.rm=TRUE)
      if(max_log_weights > -Inf){
        log_weights <- log_weights - max_log_weights
        weights <- exp(log_weights)
        ##guide_fun <- exp(log_guide_fun)
        xx <- tryCatch(
            .Call('girf_computations',
                  x=X,
                  params=params,
                  Np=Np[nt+1],
                  trackancestry=TRUE, 
                  doparRS=FALSE, 
                  weights=weights,
                  lgps=log_guide_fun,
                  fsv=array(0,dim=c(U, lookahead_steps, Np[1])), # bootgirf doesn't use fsv, set to an arbitrary val.
                  tol=tol
                  ),
            error = function (e) {
                stop(ep,conditionMessage(e),call.=FALSE) # nocov
            }
        )
        guidesim_index = guidesim_index[xx$ancestry] # update guidesim index
        cond.loglik[nt+1, s] <- xx$loglik + max_log_weights
        x <- xx$states
        log_filter_guide_fun <- xx$logfilterguides
        fcst_samp_var <- xx$newfsv
      }
      else{
        cond.loglik[nt+1, s] <- -Inf
        x <- X
        log_filter_guide_fun <- log(tol)
      }
    }
  }
  new(
    "girfd_spatPomp",
    object,
    Ninter=Ninter,
    Nguide=Nguide,
    lookahead=lookahead,
    cond.loglik=cond.loglik,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}
