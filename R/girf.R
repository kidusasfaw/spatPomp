##' Guided intermediate resampling filter (GIRF)
##'
##' An implementation of the algorithm of Park and Ionides (2019),
##' following the pseudocode in Asfaw, Ionides and King (2019).
##'
##' @name girf
##' @rdname girf
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
##' # Run GIRF
##' girfd.b <- girf(b,
##'                 Np = 100,
##'                 Ninter = length(spat_units(b)),
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
##' Upon successful completion, \code{girf} returns an object of class
##' \sQuote{girfd_spatPomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{logLik}}}{ yields an unbiased estimate of the log-likelihood of the data under the model. }
##' }
##'
##' @references
##' \park2019
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
  "girf",
  function (object, ...)
    standardGeneric("girf")
)

setMethod(
  "girf",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("girf","data")
  }
)

setMethod(
  "girf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("girf",object)
  }
)

##' @name girf-spatPomp
##' @aliases girf,spatPomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
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

    tryCatch(
      girf.internal(
        object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        params,
        tol,
        ...
      ),
      error = function (e) pomp:::pStop("girf",conditionMessage(e))
    )
  }
)

##' @name girf-girfd_spatPomp
##' @aliases girf,girfd_spatPomp-method
##' @rdname girf
##' @export
setMethod(
  "girf",
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

    girf(as(object,"spatPomp"),
         Np=Np,
         Ninter=Ninter,
         lookahead=lookahead,
         Nguide=Nguide,
         params = params,
         tol=tol,
         ...)

  }
)

girf.internal <- function (object,
        Np,
        Ninter,
        lookahead,
        Nguide,
        params,
        tol,
        ...,
        .gnsi = TRUE) {

  verbose <- FALSE
  ep <- paste0("in ",sQuote("girf"),": ")

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
  filter_guide_fun <- array(1, dim = Np[1])
  for (nt in 0:(ntimes-1)) { ## main loop
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    # get a matrix with nguides times nreps columns to propagate using rprocess
    x_with_guides <- x[,rep(1:Np[1], rep(Nguide, Np[1]))]
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
                   params=params,.gnsi=gnsi)
    xx <- tryCatch(
      .Call('do_fcst_samp_var',
            object=object,
            X=Xg,
            Np = as.integer(Np[1]),
            times=times[(nt+2):(nt+1+lookahead_steps)],
            params=params,
            gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    fcst_samp_var <- xx
    dim(fcst_samp_var) <- c(length(spat_units(object)), lookahead_steps, Np[1])

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
        skel <- pomp::flow(object, x0=X.start, t0=tt[s+1], params=params.matrix, times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],...)
        #print("skel before adjustment")
        #print(skel)
        if(s>1 && length(znames) > 0){
          skel.start <- skel[,,1]
          X.start.znames <- X.start[znames,]
          skel.start.znames <- skel.start[znames,]
          skel.end.znames <- X.start.znames + skel.start.znames
          skel[znames,,1] <- skel.end.znames
        }
        #print("skel after adjustment")
        #print(skel)
        #if(nt >= 3) return(1)
      } else {
        skel <- X
      }
      meas_var_skel <- tryCatch(
        .Call('do_theta_to_v',
              object=object,
              skel=skel,
              Np = as.integer(Np[1]),
              times=times[(nt+2):(nt+1+lookahead_steps)],
              params=params,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      dim(meas_var_skel) <- c(length(spat_units(object)), lookahead_steps, Np[1])

      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(l in 1:lookahead_steps) fcst_var_upd[,l,] <- fcst_samp_var[,l,]*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1])
      inflated_var <- meas_var_skel + fcst_var_upd
      array.params <- array(params, dim = c(length(params), length(spat_units(object)), Np[1], lookahead_steps), dimnames = list(params = names(params)))
      mmp <- tryCatch(
        .Call('do_v_to_theta',
              object=object,
              X=skel,
              vc=inflated_var,
              Np = as.integer(Np[1]),
              times=times[(nt+2):(nt+1+lookahead_steps)],
              params=array.params,
              gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      mom_match_param <- mmp
      dim(mom_match_param) <- c(length(params), length(spat_units(object)), lookahead_steps, Np[1])
      dimnames(mom_match_param) <- list(param = names(params))

      # guide functions as product (so base case is 1)
      log_guide_fun = vector(mode = "numeric", length = Np[1])

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead_steps <= 0) discount_denom_init = object@t0
        else discount_denom_init = times[nt+1+l - lookahead_steps]
        discount_factor = 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)
        # print(times[nt+1+l] - tt[s+1])
        dmeas_weights <- tryCatch(
          (vec_dmeasure(
            object,
            y=object@data[,nt+l,drop=FALSE],
            x=skel[,,l,drop = FALSE],
            times=times[nt+1+l],
            params=mom_match_param[,,l,],
            log=FALSE,
            .gnsi=gnsi
          )),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        log_dmeas_weights <- log(dmeas_weights)
        log_resamp_weights <- apply(log_dmeas_weights[,,1,drop=FALSE], 2, function(x) sum(x))*discount_factor
        log_guide_fun = log_guide_fun + log_resamp_weights
      }
      log_guide_fun[log_guide_fun < log(tol)] <- log(tol)
      log_s_not_1_weights <- log_guide_fun - log(filter_guide_fun)
      if (!(s==1 & nt!=0)){
        log_weights <- log_s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        meas_weights <- tryCatch(
          (dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=params,
            log=FALSE,
            .gnsi=gnsi
          )),
          error = function (e) {
            stop(ep,"error in calculation of dmeas_weights: ",
                 conditionMessage(e),call.=FALSE)
          }
        )
        log_meas_weights = log(meas_weights)
        gnsi <- FALSE
        log_weights <- as.numeric(log_meas_weights) + log_s_not_1_weights
      }
      max_log_weights <- max(log_weights)
      log_weights <- log_weights - max_log_weights
      weights <- exp(log_weights)
      guide_fun <- exp(log_guide_fun)
      xx <- tryCatch(
        .Call('girf_computations',
              x=X,
              params=params,
              Np=Np[nt+1],
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
      cond.loglik[nt+1, s] <- xx$loglik + max_log_weights
      # if(nt > 7 & nt < 11 & s == 1){
      # print("nt")
      # print(nt)
      # print("log guide fun")
      # print(log_guide_fun)
      # print("discount factor")
      # print(discount_factor)
      # print("filter guide fun")
      # print(filter_guide_fun)
      # print("log weights")
      # print(log_weights)
      # print("dmeas_weights")
      # print(dmeas_weights)
      # print("log_dmeas_weights")
      # print(log_dmeas_weights)
      # print("meas_weights")
      # print(meas_weights)
      # print("log_meas_weights")
      # print(log_meas_weights)
      # print("log_s_not_1_weights")
      # print(log_s_not_1_weights)
      # }
      x <- xx$states
      filter_guide_fun <- xx$filterguides
      params <- xx$params[,1]
      fcst_samp_var <- xx$newfsv
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
