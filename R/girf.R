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
##' @param h A user-provided function taking two named arguments: \code{state.vec} (representing the latent state)
##' and \code{param.vec} (representing a parameter vector for the model). It should return a scalar approximation
##' to the expected observed value given a latent state and parameter vector.
##' For more information, see the examples section below.
##' @param theta.to.v A user-provided function taking two named arguments:
##' \code{meas.mean} (representing the observation mean given a latent state - as computed using the \code{h} function above)
##' and \code{param.vec} (representing a parameter vector for the model). It should return a scalar approximation
##' to the variance of the observed value given a latent state and parameter vector.
##' For more information, see the examples section below.
##' @param v.to.theta A user-provided function taking three named arguments:
##' \code{var} (representing an empirical variance), \code{state.vec} (representing a latent state) and \code{param.vec}
##'  (representing a parameter vector for the model). The function should return a parameter vector having observation
##'   noise consistent with variance \code{var} at latent state \code{state.vec} with other parameters given by \code{param.vec}.
##' @param Ninter the number of intermediate resampling timepoints.
##' @param lookahead The number of future observations included in the guide function.
##' @param Nguide The number of simulations used to estimate state process uncertainty for each particle.
##' @param tol If the guide functions become too small (beyond floating-point precision limits), we set them to this value.
##'
##' @examples
##' # Create a simulation of a BM using default parameter set
##' b <- bm(U=3, N=10)
##'
##' # Specify the expected measurement, theta.to.v and v.to.theta
##' girf.h <- function(state.vec, param.vec){
##'   # find index matching unit_statename
##'   ix<-grep('X',names(state.vec))
##'   state.vec[ix]
##' }
##' girf.theta.to.v <- function(meas.mean, param.vec){
##'   param.vec['tau']^2
##' }
##' girf.v.to.theta <- function(var, state.vec, param.vec){
##'   param.vec['tau'] <- sqrt(var)
##'   param.vec
##' }
##' # Run GIRF
##' girfd.b <- girf(b,
##'                 Np = 100,
##'                 Ninter = length(spat_units(b)),
##'                 lookahead = 1,
##'                 Nguide = 50,
##'                 h = girf.h,
##'                 theta.to.v = girf.theta.to.v,
##'                 v.to.theta = girf.v.to.theta
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
    h = "function",
    theta.to.v = "function",
    v.to.theta = "function",
    cond.loglik="array",
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
    cond.loglik=array(data=numeric(0),dim=c(0,0)),
    Np=as.integer(NA),
    tol=as.double(NA),
    loglik=as.double(NA)
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
    h,
    theta.to.v,
    v.to.theta,
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
        h,
        theta.to.v,
        v.to.theta,
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
            h,
            theta.to.v,
            v.to.theta,
            params,
            tol,
            ...
            ) {
    if (missing(Np)) Np <- object@Np
    if (missing(tol)) tol <- object@tol
    if (missing(Ninter)) Ninter <- object@Ninter
    if (missing(Nguide)) Nguide <- object@Nguide
    if (missing(lookahead)) lookahead <- object@lookahead
    if (missing(h)) h <- object@h
    if (missing(theta.to.v)) theta.to.v <- object@theta.to.v
    if (missing(v.to.theta)) v.to.theta <- object@v.to.theta
    if (missing(params)) params <- coef(object)

    girf(as(object,"spatPomp"),
         Np=Np,
         Ninter=Ninter,
         lookahead=lookahead,
         Nguide=Nguide,
         h=h,
         theta.to.v=theta.to.v,
         v.to.theta=v.to.theta,
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
        h,
        theta.to.v,
        v.to.theta,
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
  cond.loglik <- array(0, dim = c(ntimes, Ninter))
  # initialize filter guide function
  filter_guide_fun <- array(1, dim = Np[1])
  ## begin multi-thread code
  #
  # foreach now registered outside of girf
  # doParallel::registerDoParallel(cores = NULL)
  #
  mcopts <- list(set.seed=TRUE)
  acomb <- function(...) abind::abind(..., along=3)
  acombb <- function(...) abind::abind(..., along=4)
  for (nt in 0:(ntimes-1)) { ## main loop
    # intermediate times. using seq to get S+1 points between t_n and t_{n+1} inclusive
    tt <- seq(from=times[nt+1],to=times[nt+2],length.out=Ninter+1)
    lookahead_steps = min(lookahead, ntimes-nt)
    ## for each particle get K guide particles, and fill in sample variance over K for each (lookahead value - unit - particle) combination
    fcst_samp_var <- foreach::foreach(i=1:Np[1],
         .packages=c("pomp","spatPomp"),
         .combine = 'acomb', .multicombine = TRUE,
         .options.multicore=mcopts) %dopar%  {
      Xg = array(0, dim=c(length(statenames), Nguide, lookahead_steps), dimnames = list(nvars = statenames, ng = NULL, lookahead = 1:lookahead_steps))
      xp = matrix(x[,i], nrow = nrow(x), ncol = Nguide, dimnames = list(nvars = statenames, ng = NULL))
      # get all the guides for this particle
      Xg <- rprocess(object, x0=xp, t0=times[nt+1], times=times[(nt+2):(nt+1+lookahead_steps)],
                 params=params,.gnsi=gnsi)
      fsv <- array(0, dim = c(length(spat_units(object)),lookahead_steps))
      for(u in 1:length(spat_units(object))){
        snames <- paste0(object@unit_statenames,u)
        hXgp <- apply(Xg[snames,,, drop = FALSE], MARGIN = c(2,3), FUN = h, param.vec = coef(object))
        fsv[u,] <- apply(hXgp, MARGIN = 2, FUN = var)
      }
      fsv
         }
    #print(paste0("fcst_samp_var"))
    #print(fcst_samp_var)

    # tt has S+1 (or Ninter+1) entries
    for (s in 1:Ninter){
      # get prediction simulations
      X <- rprocess(object,x0=x, t0 = tt[s], times= tt[s+1],
                    params=params,.gnsi=gnsi)
      # print(paste0("X"))
      # print(X)
      # X is now a nvars by nreps by 1 array
      X.start <- X[,,1]
      if(tt[s+1] < times[nt + 1 + lookahead_steps]){
        skel <- pomp::flow(object, x0=X.start, t0=tt[s+1], params=params.matrix, times = times[(nt + 1 + 1):(nt + 1 + lookahead_steps)],...)
      } else {
        skel <- X
      }
       #print(paste0("skel done"))
       #print(skel)
      # create measurement variance at skeleton matrix
      meas_var_skel <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))

      for(u in 1:length(object@units)){
        snames = paste0(object@unit_statenames,u)
        for (l in 1:lookahead_steps){
          hskel <- sapply(1:Np[1], function(i) apply(X=skel[snames,i,l, drop = FALSE], MARGIN = c(2,3), FUN = h, param.vec = coef(object)))
          dim(hskel) <- c(1,Np[1],1)
          meas_var_skel[u,l,] <- sapply(1:Np[1], function(i) theta.to.v(meas.mean = hskel[1,i,1], param.vec = coef(object)))
        }
      }
      #
      # CURRENT ASSUMPTIONS ARE THAT theta.to.v(state.vec, param.vec) has state.vec being a scalar
      # on the measurement scale, produced by applying h() to the skeleton.
      # This differs from the pseudocode, where state.vec would be the skeleton itself.
      #
      # h(state.vec, param.vec) should map a state vector to a scalar.
      #
      # v.to.theta(var, param.vec, state.vec) is scalar-valued.
      # it should take a scalar "var",
      # a parameter vector and a state vector -- unlike theta.to.v which wants
      # the state.vec argument to be h(state).
      #
      # there is room for code improvements here.
      #

      #
      # print(paste0("meas_var_skel"))
      # print(meas_var_skel)
      fcst_var_upd <- array(0, dim = c(length(object@units), lookahead_steps, Np[1]))
      for(u in 1:length(object@units)){
        for(l in 1:lookahead_steps){
          fcst_var_upd[u,l,] <- apply(fcst_samp_var[u,l,,drop = FALSE], MARGIN = 1,
                                      FUN = function(x) x*(times[nt+1+l] - tt[s+1])/(times[nt+1+l] - times[nt+1]))
        }
      }
      # print(paste0("fcst_var_upd"))
      # print(fcst_var_upd)
      mom_match_param <- array(0, dim = c(length(params), length(object@units), lookahead_steps, Np[1]), dimnames = list(params = names(params),unit = NULL ,lookahead = NULL, J = NULL))
      inflated_var <- meas_var_skel + fcst_var_upd
      # print(paste0("inflated_var"))
      # print(inflated_var)
      mom_match_param <- foreach::foreach(i=1:Np[1],
           .packages=c("pomp","spatPomp"),
           .combine = acombb, .multicombine = TRUE,
           .options.multicore=mcopts) %dopar%  {
        mmp <- array(0, dim = c(length(params), length(spat_units(object)), lookahead_steps), dimnames = list(params = names(params),unit = NULL ,lookahead = NULL))
        for(u in 1:length(object@units)){
          snames = paste0(object@unit_statenames,u)
          for (l in 1:lookahead_steps){
            mmp[,u,l] = v.to.theta(var = inflated_var[u,l,i], param.vec = coef(object), state.vec = skel[snames, i, l])
          }
        }
        mmp
      }
      # for(u in 1:length(object@units)){
      #   snames = paste0(object@unit_statenames,u)
      #   for (l in 1:lookahead_steps){
      #     for(p in 1:Np[1]){
      #       mom_match_param[, u, l, p] = v.to.theta(inflated_var[u,l,p], coef(object), skel[snames, p, l])
      #     }
      #   }
      # }
      # print(paste0("mom_match_param"))
      # print(mom_match_param)
      #return(1)
      # guide functions as product (so base case is 1)
      log_guide_fun = vector(mode = "numeric", length = Np[1])

      for(l in 1:lookahead_steps){
        log_dmeas_weights <- tryCatch(
          log(vec_dmeasure(
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
        log_resamp_weights <- apply(log_dmeas_weights[,,1,drop=FALSE], 2, function(x) sum(x))
        log_guide_fun = log_guide_fun + log_resamp_weights
      }
      #print("dmeas_weights")
      #print(dmeas_weights)
      log_guide_fun[log_guide_fun < log(tol)] <- log(tol)
      #print(paste0("guide_fun"))
      #print(guide_fun)
      log_s_not_1_weights <- log_guide_fun - log(filter_guide_fun)
      if (!(s==1 & nt!=0)){
        log_weights <- log_s_not_1_weights
      }
      else {
        x_3d <- x
        dim(x_3d) <- c(dim(x),1)
        rownames(x_3d)<-rownames(x)
        log_meas_weights <- tryCatch(
          log(dmeasure(
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
    h = h,
    theta.to.v = theta.to.v,
    v.to.theta = v.to.theta,
    cond.loglik=cond.loglik,
    Np=Np[1],
    tol=tol,
    loglik=sum(cond.loglik)
  )
}
