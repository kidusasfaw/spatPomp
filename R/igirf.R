##' Iterated guided intermediate resampling filter (IGIRF)
##'
##' An implementation of a parameter estimation algorithm combining
##' the intermediate resampling scheme of the guided intermediate resampling filter of Park and Ionides (2020)
##' and the parameter perturbation scheme of Ionides et al. (2015) following the pseudocode in Asfaw, et al. (2020).
##'
##' @name igirf
##' @rdname igirf
##' @include spatPomp_class.R spatPomp.R girf.R iter_filter.R
##' @family particle filter methods
##' @family spatPomp parameter estimation methods
##' @importFrom stats weighted.mean
##' @importFrom utils head
##' @inheritParams girf
##' @inheritParams pomp::mif2
##'
##' @param data an object of class \code{spatPomp} or \code{igirfd_spatPomp}
##' @param Ngirf the number of iterations of parameter-perturbed GIRF.
##'
##' @examples
##' # Complete examples are provided in the package tests
##' \dontrun{
##' igirf(bm(U=2,N=4),Ngirf=2,rw.sd = rw.sd(rho=0.02,X1_0=ivp(0.02)),cooling.type="geometric",cooling.fraction.50=0.5,Np=10,Ninter=2,lookahead=1,Nguide=5)
##' }
##' @return Upon successful completion, \code{igirf()} returns an object of class
##' \sQuote{igirfd_spatPomp}. This object contains the convergence record of the iterative algorithm with
##' respect to the likelihood and the parameters of the model (which can be accessed using the \code{traces}
##' attribute) as well as a final parameter estimate, which can be accessed using the \code{coef()}. The
##' algorithmic parameters used to run \code{igirf()} are also included.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{coef}}}{ gives the Monte Carlo maximum likelihood parameter estimate. }
##' }
##'
##' @references
##' \park2020
##'
##' \asfaw2020
NULL

rw.sd <- safecall

setClass(
  "igirfd_spatPomp",
  contains="girfd_spatPomp",
  slots=c(
    Ngirf = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix'
  )
)

setGeneric(
  "igirf",
  function (data, ...)
    standardGeneric("igirf")
)

##' @name igirf-missing
##' @aliases igirf,missing-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="missing"),
  definition=function (...) {
    stop("igirf: ","data"," is a required argument.")
  }
)

##' @name igirf-ANY
##' @aliases igirf,ANY-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop("igirf is undefined for ", sQuote(data), "of class ", sQuote(class(data)), ".")
  }
)

##' @name igirf-spatPomp
##' @aliases igirf,spatPomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="spatPomp"),
  definition=function (data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
    Ninter,lookahead=1,Nguide,kind=c('bootstrap', 'moment'),
    tol = 1e-300,
    ..., verbose = getOption("verbose", FALSE)) {
    if(missing(Ninter)) Ninter <- length(unit_names(data))
    kind = match.arg(kind)
    tryCatch(
      igirf.internal(data,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
        Ninter,lookahead,Nguide,kind,tol = tol,
        ...,verbose=verbose),
      error = function (e) stop("igirf: ", e)
    )
  }
)

##' @name igirf-igirfd_spatPomp
##' @aliases igirf,igirfd_spatPomp-method
##' @rdname igirf
##' @export
setMethod(
  "igirf",
  signature=signature(data="igirfd_spatPomp"),
  function (data,Ngirf,Np,rw.sd,cooling.type, cooling.fraction.50, Ninter,
    lookahead,Nguide,kind=c('bootstrap','moment'),tol, ...,
    verbose = getOption("verbose", FALSE)) {
    if (missing(Ngirf)) Ngirf <- data@Ngirf
    if (missing(rw.sd)) rw.sd <- data@rw.sd
    if (missing(cooling.type)) cooling.type <- data@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- data@cooling.fraction.50
    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol
    if (missing(Ninter)) Ninter <- data@Ninter
    if (missing(Nguide)) Nguide <- data@Nguide
    if (missing(lookahead)) lookahead <- data@lookahead

    igirf(as(data,"spatPomp"), Ngirf=Ngirf, Np=Np,rw.sd = rw.sd, cooling.type = cooling.type,
      cooling.fraction.50 = cooling.fraction.50, tol=tol, Ninter=Ninter, Nguide=Nguide,
      kind=kind, lookahead=lookahead, ..., verbose=verbose)
  }
)

igirf.internal <- function (object,Ngirf,Np,rw.sd,cooling.type,cooling.fraction.50,
  Ninter,lookahead,Nguide,kind,
  tol, ...,
  .ndone = 0L, .indices = integer(0),.paramMatrix = NULL,.gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)
  p_object <- pomp(object,...,verbose=FALSE)
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
    unit_obsnames = object@unit_obsnames)
  if (undefined(object@rprocess) || undefined(object@dmeasure))
    stop(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Ngirf) != 1 || !is.numeric(Ngirf) || !is.finite(Ngirf) || Ngirf < 1)
    stop(sQuote("Ngirf")," must be a positive integer.")
  Ngirf <- as.integer(Ngirf)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))

  if (is.null(Np)) {
    stop(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        stop("if ",sQuote("Np"),
          " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    stop(sQuote("Np"),
      " must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1) {
    Np <- rep(Np,times=ntimes)
  } else if (length(Np) > ntimes) {
    if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
      warning("in igirf ","Np[k] ignored for k > ",sQuote("length(time(object))"),".")
    }
    Np <- head(Np,ntimes)
  } else if (length(Np) < ntimes) {
    stop(sQuote("Np")," must have length 1 or ",
      sQuote("length(time(object))"),".")
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    stop(sQuote("Np")," must be a positive integer.")

  Np <- as.integer(Np)
  Np <- c(Np,Np[1L])

  if (missing(rw.sd))
    stop(sQuote("rw.sd")," must be specified!")
  ##rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
  rw.sd <- perturbn.kernel.sd(rw.sd,
    time=igirf_rw_sd_times(times=time(object),Ninter=Ninter),
    paramnames=names(start))
  if (missing(cooling.fraction.50))
    stop(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
        !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
          cooling.fraction.50 > 1)
    stop(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
      dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }

  traces <- array(dim=c(Ngirf+1,length(start)+1),
    dimnames=list(iteration=seq.int(.ndone,.ndone+Ngirf),
      variable=c('loglik',names(start))))
  traces[1L,] <- c(NA,start)

  pompLoad(object,verbose=FALSE)
  on.exit(pompUnload(object,verbose=FALSE))

  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
    .gnsi=gnsi)

  ## iterate the filtering
  for (n in seq_len(Ngirf)) {
    if(kind == 'moment'){
      g <- igirf.momgirf(object=object,Ninter=Ninter,Nguide=Nguide,lookahead=lookahead,
        params=paramMatrix,
        Np=Np,girfiter=.ndone+n,cooling.fn=cooling.fn,rw.sd=rw.sd,tol=tol,
        verbose=verbose,.indices=.indices, .gnsi=gnsi)
    }
    if(kind == 'bootstrap'){
      g <- igirf.bootgirf(
        object=object,Ninter=Ninter,Nguide=Nguide,lookahead=lookahead,
        params=paramMatrix,
        Np=Np,girfiter=.ndone+n,cooling.fn=cooling.fn,rw.sd=rw.sd,tol=tol,
        verbose=verbose,.indices=.indices, .gnsi=gnsi
      )
    }
    gnsi <- FALSE

    paramMatrix <- g@paramMatrix
    traces[n+1,-1] <- coef(g)
    traces[n+1,1] <- g@loglik
    .indices <- .indices

    if (verbose) {
      cat("igirf iteration",n,"of",Ngirf,"completed with likelihood ", g@loglik, "\n")
      print(coef(g))
    }
  }

  g@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
    .gnsi=gnsi)

  new(
    "igirfd_spatPomp",
    g,
    Ngirf=Ngirf,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )
}

igirf.momgirf <- function (object, params, Ninter, lookahead, Nguide,
  Np, girfiter, rw.sd, cooling.fn, tol,
  ..., verbose, .indices = integer(0), .gnsi = TRUE) {

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  girfiter <- as.integer(girfiter)
  Np <- as.integer(Np)
  ep <- paste0("in ",sQuote("igirf"),": ")

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop(sQuote("tol")," should be a small positive number.")

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    stop(sQuote(".indices")," has improper length.")

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))

  cond.loglik <- array(0, dim = c(ntimes, Ninter))

  znames <- object@accumvars

                                        # Initialize filter guide function
  log_filter_guide_fun <- array(0, dim = Np[1])
  for (nt in 0:(ntimes-1)) {
    pmag <- cooling.fn(nt,girfiter)$alpha*rw.sd[,nt]
    params <- .Call(randwalk_perturbation_spatPomp,params,pmag)
    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
    if (nt == 0) {
      x <- rinit(object,params=tparams)
      statenames <- rownames(x)
    }
    ## Intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps <- min(lookahead, ntimes-nt)
    ## Four-dimensional array: nvars by nguide by ntimes by nreps
    Xg <- array(0,
      dim=c(length(statenames), Nguide, lookahead_steps, Np[1]),
      dimnames = list(
        nvars = statenames,
        ng = NULL,
        lookahead = seq_len(lookahead_steps),
        nreps = NULL
      )
    )
                                        # For each particle get K guide particles, and fill in sample variance over K for each (lookahead value - unit - particle) combination
    fcst_samp_var <- array(0, dim = c(length(unit_names(object)), lookahead_steps, Np[1]))
    x_with_guides <- x[,rep(1:Np[1], rep(Nguide, Np[1]))]
    tp_with_guides <- tparams[,rep(1:Np[1], rep(Nguide, Np[1]))]
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
      params=tp_with_guides,.gnsi=gnsi)
    xx <- tryCatch(
      .Call(do_fcst_samp_var,
        object=object,
        X=Xg,
        Np = as.integer(Np[1]),
        times=times[(nt+2):(nt+1+lookahead_steps)],
        params=tp_with_guides,
        gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    fcst_samp_var <- xx
    dim(fcst_samp_var) <- c(length(unit_names(object)), lookahead_steps, Np[1])

    for (s in 1:Ninter){
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
                                        # Get prediction simulations: nvars by nreps by 1 array
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
        params=tparams,.gnsi=gnsi)
      if(s>1 && length(znames)>0){
        x.znames <- x[znames,]; dim(x.znames) <- c(dim(x.znames),1)
        X[znames,,] <- X[znames,,,drop=FALSE] + x.znames
      }
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp::flow(object,
          x0=X.start,
          t0=tt[s+1],
          params=tparams,
          times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],
          method = 'adams')
        skel.start <- skel[,,1]
        X.start.znames <- X.start[znames,]
        skel.start.znames <- skel.start[znames,]
        skel.end.znames <- X.start.znames + skel.start.znames
        skel[znames,,1] <- skel.end.znames
      } else {
        skel <- X
      }
      meas_var_skel <- tryCatch(
        .Call(do_theta_to_v,
          object=object,
          X=skel,
          Np = as.integer(Np[1]),
          times=times[(nt+2):(nt+1+lookahead_steps)],
          params=tparams,
          gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      dim(meas_var_skel) <- c(length(unit_names(object)), lookahead_steps, Np[1])
      fcst_var_upd <- array(0, dim = c(length(unit_names(object)), lookahead_steps, Np[1]))
      for(l in 1:lookahead_steps) fcst_var_upd[,l,] <- fcst_samp_var[,l,]*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1])
      inflated_var <- meas_var_skel + fcst_var_upd
      array.tparams <- array(NA, dim = c(dim(tparams)[1], length(unit_names(object)), Np[1], lookahead_steps), dimnames = list(tparams = rownames(tparams)))
      for(i in 1:length(unit_names(object))) array.tparams[,i,,1:lookahead_steps] <- tparams
      mmp <- tryCatch(
        .Call(do_v_to_theta,
          object=object,
          X=skel,
          vc=inflated_var,
          Np = as.integer(Np[1]),
          times=times[(nt+2):(nt+1+lookahead_steps)],
          params=array.tparams,
          gnsi=TRUE),
        error = function (e) {
          stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
      )
      mom_match_param <- mmp
      dim(mom_match_param) <- c(dim(tparams)[1], length(unit_names(object)), lookahead_steps, Np[1])
      dimnames(mom_match_param) <- list(tparam = rownames(tparams))
      log_guide_fun <- vector(mode = "numeric", length = Np[1]) + 1

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead <= 0) discount_denom_init <- object@t0
        else discount_denom_init <- times[nt+1+l - lookahead]
        discount_factor <- 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)/ifelse(lookahead==1,2,1) ## to ensure that the discount factor does not become too small for L=1 and small s (which can lead to very uninformative guide function), increase the discount factor to at least 1/2 when L=1.
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
        if(any(is.na(log_dmeas_weights))){
                                        # find particle with the NA
          na_ix <- which(is.na(log_dmeas_weights[,,1]))[1]
          na_ix_col <- (na_ix %/% nunits) + (na_ix %% nunits > 0)
          illegal_dunit_measure_error(
            time=times[nt+1+l],
            lik=log_dmeas_weights[,na_ix_col,1,drop=FALSE],
            datvals=object@data[,nt+l],
            states=skel[,na_ix_col,l],
            params=mom_match_param[,,l,na_ix_col]
          )
        }
        log_resamp_weights <- apply(log_dmeas_weights[,,1,drop=FALSE], 2, function(x) sum(x))*discount_factor
        log_guide_fun <- log_guide_fun + log_resamp_weights
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
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=tparams,
            log=TRUE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of log_meas_weights: ",
              conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        log_weights <- as.numeric(log_meas_weights) + log_s_not_1_weights
      }
      if (nt == ntimes-1 & s==Ninter) {
        if (any(log_weights>-Inf)) {
          coef(object,transform=TRUE) <- apply(params,1L,weighted.mean,w=exp(log_weights))
        } else {
          warning("igirf: ","filtering failure at last filter iteration; using ",
            "unweighted mean for point estimate.")
          coef(object,transform=TRUE) <- apply(params,1L,mean)
        }
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
            doparRS=TRUE,
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
        params <- xx$params
        fcst_samp_var <- xx$newfsv
      }
      else{
        cond.loglik[nt+1, s] <- log(tol)
        x <- X[,,1]
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
    paramMatrix=params,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}

igirf.bootgirf <- function (object, params, Ninter, lookahead, Nguide,
  Np, girfiter, rw.sd, cooling.fn, tol,
  ..., verbose, .indices = integer(0), .gnsi = TRUE) {

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  girfiter <- as.integer(girfiter)
  Np <- as.integer(Np)
  ep <- paste0("in ",sQuote("ibootgirf"),": ")

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop(ep,sQuote("tol")," should be a small positive number.",call.=FALSE)

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    stop(ep,sQuote(".indices")," has improper length.",call.=FALSE)

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  nunits <- length(unit_names(object))

  cond.loglik <- array(0, dim = c(ntimes, Ninter))

  znames <- object@accumvars

  log_filter_guide_fun <- array(0, dim = Np[1])
  for (nt in 0:(ntimes-1)) {
    pmag <- cooling.fn(nt,girfiter)$alpha*rw.sd[,nt]
    params <- .Call(randwalk_perturbation_spatPomp,params,pmag)
    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
    if (nt == 0) {
      x <- rinit(object,params=tparams)
      statenames <- rownames(x)
    }
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps <- min(lookahead, ntimes-nt)
    ## Four-dimensional array: nvars by nguide by ntimes by nreps
    Xg <- array(0,
      dim=c(length(statenames), Nguide, lookahead_steps, Np[1]),
      dimnames = list(
        nvars = statenames,
        ng = NULL,
        lookahead = seq_len(lookahead_steps),
        nreps = NULL
      )
    )
    guidesim_index <- seq_len(Np[1]) ## the index for guide simulations (to be updated each time resampling occurs)
    x_with_guides <- x[,rep(guidesim_index, each = Nguide)]
    tp_with_guides <- tparams[,rep(guidesim_index, each = Nguide)]
    Xg <- rprocess(object, x0=x_with_guides, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
      params=tp_with_guides,.gnsi=gnsi)
    Xskel <- pomp::flow(object,
      x0=x,
      t0=times[nt+1],
      params=tparams,
      times = times[(nt+2):(nt+1+lookahead_steps)],
      method = 'adams')
    resids <- Xg - Xskel[,rep(1:Np[1], each=Nguide),,drop=FALSE] # residuals
    rm(Xg, Xskel, x_with_guides)
    for (s in 1:Ninter){
      tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)
      tp_with_guides <- tparams[,rep(1:Np[1], rep(Nguide, Np[1]))]
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
        params=tparams,.gnsi=gnsi)
      if(s>1 && length(znames)>0){
        x.znames <- x[znames,]; dim(x.znames) <- c(dim(x.znames),1)
        X[znames,,] <- X[znames,,,drop=FALSE] + x.znames
      }
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp::flow(object,
          x0=X.start,
          t0=tt[s+1],
          params=tparams,
          times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],
          method = 'adams')
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
      log_guide_fun <- vector(mode = "numeric", length = Np[1]) + 1

      for(l in 1:lookahead_steps){
        if(nt+1+l-lookahead <= 0) discount_denom_init <- object@t0
        else discount_denom_init <- times[nt+1+l - lookahead]
        discount_factor <- 1 - (times[nt+1+l] - tt[s+1])/(times[nt+1+l] - discount_denom_init)/ifelse(lookahead==1,2,1) ## to ensure that the discount factor does not become too small for L=1 and small s (which can lead to very uninformative guide function), increase the discount factor to at least 1/2 when L=1.

                                        # construct pseudo-simulations by adding simulated noise terms (residuals) to the skeletons
        pseudosims <- skel[,rep(1:Np[1], each=Nguide),l,drop=FALSE] +
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),l,drop=FALSE] -
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),1,drop=FALSE] +
          resids[,rep(guidesim_index-1, each=Nguide)*Nguide+rep(1:Nguide, Np[1]),1,drop=FALSE] * sqrt((times[nt+2]-tt[s+1])/(times[nt+2]-times[nt+1]))
        rm(skel)
        log_dmeas_weights <- tryCatch(
        (vec_dmeasure(
          object,
          y=object@data[,nt+l,drop=FALSE],
          x=pseudosims,
          times=times[nt+1+l],
          params=tp_with_guides,
          log=TRUE,
          .gnsi=gnsi
        )),
        error = function (e) {
          stop(ep,"error in calculation of log_dmeas_weights: ",
            conditionMessage(e),call.=FALSE)
        }
        )
        if(any(is.na(log_dmeas_weights))){
                                        # find particle with the NA
          na_ix <- which(is.na(log_dmeas_weights[,,1]))[1]
          na_ix_col <- (na_ix %/% nunits) + (na_ix %% nunits > 0)
          illegal_dunit_measure_error(
            time=times[nt+1+l],
            lik=log_dmeas_weights[,na_ix_col,1,drop=FALSE],
            datvals=object@data[,nt+l],
            states=pseudosims[,na_ix_col,1L],
            params=tp_with_guides[,na_ix_col]
          )
        }
        ldw <- array(log_dmeas_weights, c(nunits,Nguide,Np[1])) # log_dmeas_weights is an array with dim U*(Np*Nguide)*1. Reorder it as U*Nguide*Np
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
          dmeasure(
            object,
            y=object@data[,nt,drop=FALSE],
            x=x_3d,
            times=times[nt+1],
            params=tparams,
            log=TRUE,
            .gnsi=gnsi
          ),
          error = function (e) {
            stop(ep,"error in calculation of log_meas_weights: ",
              conditionMessage(e),call.=FALSE)
          }
        )
        gnsi <- FALSE
        log_weights <- as.numeric(log_meas_weights) + log_s_not_1_weights
      }
      if (nt == ntimes-1 & s==Ninter) {
        if (any(log_weights>-Inf)) {
          coef(object,transform=TRUE) <- apply(params,1L,weighted.mean,w=exp(log_weights))
        } else {
          pWarn("ibootgirf","filtering failure at last filter iteration; using ",
            "unweighted mean for point estimate.")
          coef(object,transform=TRUE) <- apply(params,1L,mean)
        }
      }
      max_log_weights <- max(log_weights, na.rm = TRUE)
      if(max_log_weights > -Inf){
        log_weights <- log_weights - max_log_weights
        weights <- exp(log_weights)
        xx <- tryCatch(
          .Call(girf_computations,
            x=X,
            params=params,
            Np=Np[nt+1],
            trackancestry=TRUE,
            doparRS=TRUE,
            weights=weights,
            lgps=log_guide_fun,
            fsv=array(0,dim=c(length(unit_names(object)), lookahead_steps, Np[1])), # bootgirf2 doesn't use fsv, set to an arbitrary val.
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
        params <- xx$params
      }
      else{
        cond.loglik[nt+1, s] <- log(tol)
        x <- X[,,1]
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
    paramMatrix=params,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}

illegal_dunit_measure_error <- function(time, lik, datvals, states, params){
  showvals <- c(time=time,lik=lik,datvals,states,params)
  m1 <- formatC(names(showvals),preserve.width="common")
  m2 <- formatC(showvals,digits=6,width=12,format="g",preserve.width="common")
  stop(
    sQuote("dunit_measure")," returns illegal value.\n",
    "Likelihood, data, states, and parameters are:\n",
    paste0(m1,": ",m2,collapse="\n")
  )
}

igirf_rw_sd_times <- function(times, Ninter){
  rw_sd_times <- c()
  for(i in seq(length(times)-1)) rw_sd_times <- c(
                                   rw_sd_times, seq(from=times[i],to=times[i+1],length.out=Ninter+1))
  return(rw_sd_times)
}




