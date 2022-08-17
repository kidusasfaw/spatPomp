##' Guided intermediate resampling filter (GIRF)
##'
##' An implementation of the algorithm of Park and Ionides (2020),
##' following the pseudocode in Asfaw et al. (2020).
##'
##' @name girf
##' @rdname girf
##' @include spatPomp_class.R spatPomp.R
##' @author Kidus Asfaw
##' @family particle filter methods
##'
##'
##' @inheritParams abf
##' @param Ninter the number of intermediate resampling time points. By default, this is set equal to the number of units.
##' @param lookahead The number of future observations included in the guide function.
##' @param Nguide The number of simulations used to estimate state process uncertainty for each particle.
##' @param kind One of two types of guide function construction. Defaults to \code{'bootstrap'}. See Park and Ionides (2020) for more details.
##' @param tol If all of the guide function evaluations become too small (beyond floating-point precision limits), we set them to this value.
##'
##' @examples
##' # Complete examples are provided in the package tests
##' \dontrun{
##' #
##' # Create a simulation of a Brownian motion
##' b <- bm(U=2, N=5)
##'
##' # Run GIRF
##' girfd_bm <- girf(b,
##'                  Np = 10,
##'                  Ninter = length(unit_names(b)),
##'                  lookahead = 1,
##'                  Nguide = 10
##' )
##' # Get the likelihood estimate from GIRF
##' logLik(girfd_bm)
##'
##' # Compare with the likelihood estimate from particle filter
##' pfd_bm <- pfilter(b, Np = 10)
##' logLik(pfd_bm)
##' }
##' @return Upon successful completion, \code{girf()} returns an object of class
##' \sQuote{girfd_spatPomp} which contains the algorithmic parameters that were used to
##' run \code{girf()} and the resulting log likelihood estimate.
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
##' \asfaw2020
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

##' @name girf-missing
##' @aliases girf,missing-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("girf: ","data"," is a required argument.")
  }
)

##' @name girf-ANY
##' @aliases girf,ANY-method
##' @rdname girf
##' @export
setMethod(
  "girf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop("girf is undefined for ", sQuote(object), "of class ", sQuote(class(object)), ".")
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
    lookahead=1,
    Nguide,
    kind=c('bootstrap','moment'),
    tol,
    ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(tol)) tol <- 1e-300
    if (missing(Ninter)) Ninter <- length(unit_names(object))
    kind = match.arg(kind)

    if(kind == 'moment'){
      tryCatch(
        g <- momgirf.internal(
          object,
          Np=Np,
          Ninter=Ninter,
          lookahead=lookahead,
          Nguide=Nguide,
          tol=tol,
          ...,
          verbose=verbose
        ),
        error = function (e) pStop("girf",conditionMessage(e))
      )
      return(g)
    }

    if(kind == 'bootstrap'){
      tryCatch(
        g <- bootgirf.internal(
          object,
          Np=Np,
          Ninter=Ninter,
          lookahead=lookahead,
          Nguide=Nguide,
          tol=tol,
          ...,
          verbose=verbose
        ),
        error = function (e) pStop("girf",conditionMessage(e))
      )
      return(g)
    }

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
    kind = c('bootstrap','moment'),
    tol,
    ...
  ) {
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    if (missing(Ninter)) Ninter <- object@Ninter
    if (missing(Nguide)) Nguide <- object@Nguide
    if (missing(lookahead)) lookahead <- object@lookahead

    girf(as(object,"spatPomp"),
      Np=Np,
      Ninter=Ninter,
      lookahead=lookahead,
      Nguide=Nguide,
      kind=kind,
      tol=tol,
      ...)

  }
)

momgirf.internal <- function (object,
  Np,
  Ninter,
  lookahead,
  Nguide,
  tol,
  ...,
  verbose,
  .gnsi = TRUE) {
  verbose <- as.logical(verbose)
  p_object <- pomp(object,...,verbose=verbose)
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

  ep <- paste0("in ",sQuote("girf"),": ")
  params <- coef(object)

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)

  params.matrix <- matrix(params,nrow=length(params), ncol = Np[1])
  rownames(params.matrix) <- names(params)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  if (missing(Np) || is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small positive number.")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  x <- init.x
  znames <- object@accumvars
  cond.loglik <- array(0, dim = c(ntimes, Ninter))
  log_filter_guide_fun <- array(0, dim = Np[1])
  for (nt in 0:(ntimes-1)) { ## main loop
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps <- min(lookahead, ntimes-nt)
                                        # Get a matrix with nguides times nreps columns to propagate using rprocess
    x_with_guides <- x[,rep(1:Np[1], rep(Nguide, Np[1]))]
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
      params=params,.gnsi=gnsi)
    xx <- tryCatch(
      .Call(do_fcst_samp_var,
        object=object,
        X=Xg,
        Np = as.integer(Np[1]),
        times=times[(nt+2):(nt+1+lookahead_steps)],
        params=params,
        gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
    fcst_samp_var <- xx
    dim(fcst_samp_var) <- c(length(unit_names(object)), lookahead_steps, Np[1])
    for (s in 1:Ninter){
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
        params=params,.gnsi=gnsi)
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
      meas_var_skel <- tryCatch(
        .Call(do_theta_to_v,
          object=object,
          X=skel,
          Np = as.integer(Np[1]),
          times=times[(nt+2):(nt+1+lookahead_steps)],
          params=params,
          gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      dim(meas_var_skel) <- c(length(unit_names(object)), lookahead_steps, Np[1])

      fcst_var_upd <- array(0, dim = c(length(unit_names(object)), lookahead_steps, Np[1]))
      for(l in 1:lookahead_steps) fcst_var_upd[,l,] <- fcst_samp_var[,l,]*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1])
      inflated_var <- meas_var_skel + fcst_var_upd
      array.params <- array(params, dim = c(length(params), length(unit_names(object)), Np[1], lookahead_steps), dimnames = list(params = names(params)))
      mmp <- tryCatch(
        .Call(do_v_to_theta,
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
      dim(mom_match_param) <- c(length(params), length(unit_names(object)), lookahead_steps, Np[1])
      dimnames(mom_match_param) <- list(param = names(params))
      log_guide_fun = vector(mode = "numeric", length = Np[1])

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead <= 0) discount_denom_init = object@t0
        else discount_denom_init = times[nt+1+l - lookahead]
        discount_factor = 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)/ifelse(lookahead==1,2,1) ## to ensure that the discount factor does not become too small for L=1 and small s (which can lead to very uninformative guide function), increase the discount factor to at least 1/2 when L=1.
        log_dmeas_weights <- tryCatch(
        (vec_dmeasure(
          object,
          y=object@data[,nt+l,drop=FALSE],
          x=skel[,,l,drop = FALSE],
          times=times[nt+1+l],
          params=mom_match_param[,,l,],
          log=TRUE,
          .gnsi=gnsi
        )),
        error = function (e) {
          stop(ep,"error in calculation of log_dmeas_weights: ",
            conditionMessage(e),call.=FALSE)
        }
        )
        log_resamp_weights <- apply(log_dmeas_weights[,,1,drop=FALSE], 2, sum)*discount_factor
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
      max_log_weights <- max(log_weights)
      if(max_log_weights > -Inf){
        log_weights <- log_weights - max_log_weights
        weights <- exp(log_weights)
        xx <- tryCatch(
          .Call(girf_computations,
            x=X,
            params=params,
            Np=Np[nt+1],
            trackancestry=FALSE,
            doparRS=FALSE,
            weights=weights,
            lgps=log_guide_fun,
            fsv=fcst_samp_var,
            tol=tol
          ),
          error = function (e) {
            stop(ep,conditionMessage(e),call.=FALSE) # nocov
          }
        )
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

bootgirf.internal <- function (object,
  Np,
  Ninter,
  lookahead,
  Nguide,
  tol,
  ...,
  verbose,
  .gnsi = TRUE) {

  ep <- paste0("in ",sQuote("girf"),": ")

  p_object <- pomp(object,...,verbose=verbose)
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

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  if (length(params)==0)
    stop(ep,sQuote("params")," must be specified",call.=FALSE)

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)

  params.matrix <- matrix(params,nrow=length(params), ncol = Np[1])
  rownames(params.matrix) <- names(params)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  U <- length(unit_names(object))
  if (missing(Np) || is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np")," is a function, it must return ",
          "a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np")," must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1)
    Np <- rep(Np,times=ntimes+1)
  else if (length(Np) != (ntimes+1))
    pStop_(sQuote("Np")," must have length 1 or length ",ntimes+1,".")

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_("number of particles, ",sQuote("Np"),", must be a positive integer.")

  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small positive number.")

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  init.x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)
  x <- init.x
  znames <- object@accumvars
  cond.loglik <- array(0, dim = c(ntimes, Ninter))
  log_filter_guide_fun <- array(0, dim = Np[1])
  for (nt in 0:(ntimes-1)) { ## main loop
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
                                        # Get a matrix with nguides times nreps columns to propagate using rprocess
    x_with_guides <- x[,rep(1:Np[1], each=Nguide)]
    guidesim_index <- 1:Np[1] # the index for guide simulations (to be updated each time resampling occurs)
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)], params=params,.gnsi=gnsi)
    Xskel <- tryCatch( # skeleton
      pomp::flow(object,
        x0=x,
        t0=times[nt+1],
        params=params.matrix,
        times = times[(nt+2):(nt+1+lookahead_steps)],
        ...),
      error = function (e) {
        pomp::flow(object,
          x0=x,
          t0=times[nt+1],
          params=params.matrix,
          times = times[(nt+2):(nt+1+lookahead_steps)],
          method = 'adams')
      }
    )
    resids <- Xg - Xskel[,rep(1:Np[1], each=Nguide),,drop=FALSE] # residuals
    rm(Xg, Xskel, x_with_guides)
    for (s in 1:Ninter){
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
        params=params,.gnsi=gnsi)
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

      log_guide_fun = vector(mode = "numeric", length = Np[1])

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead <= 0) discount_denom_init = object@t0
        else discount_denom_init = times[nt+1+l - lookahead]
        discount_factor = 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)/ifelse(lookahead==1,2,1) ## to ensure that the discount factor does not become too small for L=1 and small s (which can lead to very uninformative guide function), increase the discount factor to at least 1/2 when L=1.

                                        # construct pseudo-simulations by adding simulated noise terms (residuals) to the skeletons
        pseudosims <- skel[,rep(1:Np[1], each=Nguide),l,drop=FALSE] +
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),l,drop=FALSE] -
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),1,drop=FALSE] +
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),1,drop=FALSE] * sqrt((times[nt+2]-tt[s+1])/(times[nt+2]-times[nt+1]))

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
          .Call(girf_computations,
            x=X,
            params=params,
            Np=Np[nt+1],
            trackancestry=TRUE,
            doparRS=FALSE,
            weights=weights,
            lgps=log_guide_fun,
            fsv=array(0,dim=c(length(unit_names(object)), lookahead_steps, Np[1])), # bootgirf2 doesn't use fsv, set to an arbitrary val.
            tol=tol
          ),
          error = function (e) {
            stop(ep,conditionMessage(e),call.=FALSE) # nocov
          }
        )
        guidesim_index <- guidesim_index[xx$ancestry] # update guidesim index
        cond.loglik[nt+1, s] <- xx$loglik + max_log_weights
        x <- xx$states
        log_filter_guide_fun <- xx$logfilterguides
      }
      else{
        cond.loglik[nt+1, s] <- -Inf
        x <- X
        log_filter_guide_fun <- log(tol)
      }
      gc()
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

