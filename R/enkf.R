##' Generalized Ensemble Kalman filter (EnKF)
##'
##' A function to perform filtering using the ensemble Kalman filter of Evensen, G. (1994).
##' This function is generalized to allow for an measurement covariance matrix that varies over time.
##' This is useful if the measurement model varies with the state.
##'
##' @name enkf
##' @rdname enkf
##' @include spatPomp_class.R spatPomp.R
##' @aliases enkf  enkf,ANY-method enkf,missing-method
##' @family particle filter methods
##' @family \pkg{spatPomp} filtering methods
##' @importFrom stats rnorm
##' @inheritParams abf
##'
##' @param data A \code{spatPomp} object.
##' @param Np The number of Monte Carlo particles used to approximate the filter distribution.
##' @examples
##' # Create a simulation of a Brownian motion
##' b <- bm(U=6, N=10)
##'
##' # Run EnKF
##' enkfd_bm <- enkf(b, Np = 100)
##'
##' # Get a likelihood estimate
##' logLik(enkfd_bm)

##' @return
##' An object of class \sQuote{enkfd_spatPomp}.
##'
##' @references \Evensen1994
##'
##' \Evensen2009
##'
##' \Anderson2001
##'
NULL

setClass(
  "enkfd_spatPomp",
  contains="kalmand_pomp",
  slots=c(
    paramMatrix = 'array',
    indices = 'vector',
    unit_names = 'character',
    unit_statenames = 'character',
    unit_obsnames = 'character',
    dunit_measure = 'pomp_fun',
    runit_measure = 'pomp_fun',
    eunit_measure = 'pomp_fun',
    vunit_measure = 'pomp_fun',
    munit_measure = 'pomp_fun'
  ),
  prototype=prototype(
    paramMatrix=array(data=numeric(0),dim=c(0,0)),
    indices=integer(0),
    dunit_measure = pomp:::pomp_fun(slotname="dunit_measure"),
    runit_measure = pomp:::pomp_fun(slotname="runit_measure"),
    eunit_measure = pomp:::pomp_fun(slotname="eunit_measure"),
    vunit_measure = pomp:::pomp_fun(slotname="vunit_measure"),
    munit_measure = pomp:::pomp_fun(slotname="munit_measure")
  )
)

setMethod(
  "enkf",
  signature=signature(data="missing"),
  definition=function (...) {
    stop("enkf: ","data"," is a required argument.")
  }
)

setMethod(
  "enkf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop("enkf is undefined for ", sQuote(data), "of class ", sQuote(class(data)), ".")
  }
)

## GENERALIZED ENSEMBLE KALMAN FILTER (enkf)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$
## Forecast: $Y_t=h(X_t)$
## Forecast mean: $N_t=\langle Y \rangle$.
## Forecast variance: $S_t=\langle\langle Y \rangle\rangle$
## State/forecast covariance: $W_t=\langle\langle X,Y\rangle\rangle$
## Kalman gain: $K_t = W_t\,S_t^{-1}$
## New observation: $y_t\in \mathbb{R}^{n\times 1}$
## Updated ensemble: $X^u_{t}=X_t + K_t\,(y_t - Y_t)$
## Filter mean: $m_t=\langle X^u_t \rangle = \frac{1}{q} \sum\limits_{i=1}^q x^{u_i}_t$

##' @name enkf-spatPomp
##' @aliases enkf,spatPomp-method
##' @rdname enkf
##' @export
setMethod(
  "enkf",
  signature=signature(data="spatPomp"),
  function (data,
            Np,
            ..., verbose = getOption("verbose", FALSE)) {
    tryCatch(
      sp <- enkf.internal(
        data,
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) stop("enkf",conditionMessage(e))
    )
  }
)

enkf.internal <- function (object,
                           Np,
                           ..., verbose) {

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

  if (pomp:::undefined(object@rprocess))
    pomp:::pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")
  if (pomp:::undefined(object@eunit_measure))
    pomp:::pStop_(paste(sQuote(c("eunit_measure")),collapse=", ")," is a needed basic component.")
  if (pomp:::undefined(object@vunit_measure))
    pomp:::pStop_(paste(sQuote(c("vunit_measure")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  params <- coef(object)

  t <- time(object)
  tt <- time(object,t0=TRUE)
  ntimes <- length(t)

  y <- obs(object)
  nobs <- nrow(y)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL))
  X <- rinit(object,params=params,nsim=Np)
  nvar <- nrow(X)


  filterMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
  predMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
  forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
  condlogLik <- numeric(ntimes)

  for (k in seq_len(ntimes)) {
    ## advance ensemble according to state process
    X <- rprocess(object,x0=X,t0=tt[k],times=tt[k+1],params=params)

    # data
    yk <- y[,k]
    # ensemble of forecasts
    Y <- tryCatch(
      .Call('do_theta_to_e',
            object=object,
            X=X,
            Np = as.integer(Np[1]),
            times=tt[k+1],
            params=params,
            gnsi=TRUE),
      error = function (e) {
        stop("ep",conditionMessage(e),call.=FALSE) # nocov
      }
    )
    # variance of artificial noise (i.e. R) computed using vmeasure
    meas_var <- tryCatch(
      .Call('do_theta_to_v',
            object=object,
            X=X,
            Np = as.integer(Np[1]),
            times=tt[k+1],
            params=params,
            gnsi=TRUE),
      error = function (e) {
        stop("ep",conditionMessage(e),call.=FALSE) # nocov
      }
    )
    dim(meas_var) <- c(length(unit_names(object)),  Np[1])
    R <- diag(rowMeans(meas_var))
    sqrtR <- tryCatch(
      t(chol(R)),                     # t(sqrtR)%*%sqrtR == R
      error = function (e) {
        pomp:::pStop_("degenerate ",sQuote("R"), "at time ", sQuote(k), ": ",conditionMessage(e))
      }
    )
    X <- X[,,1]
    predMeans[,k] <- pm <- rowMeans(X) # prediction mean
    dim(Y) <- c(length(unit_names(object)), Np[1])

    # forecast mean
    ym <- rowMeans(Y)
    X <- X-pm
    Y <- Y-ym

    fv <- tcrossprod(Y)/(Np-1)+R    # forecast variance
    vyx <- tcrossprod(Y,X)/(Np-1)   # forecast/state covariance

    svdS <- svd(fv,nv=0)            # singular value decomposition
    Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d) # transpose of Kalman gain
    Ek <- sqrtR%*%matrix(rnorm(n=nobs*Np),nobs,Np) # artificial noise
    resid <- y[,k]-ym

    X <- X+pm+crossprod(Kt,resid-Y+Ek)

    condlogLik[k] <- sum(dnorm(x=crossprod(svdS$u,resid),mean=0,sd=sqrt(svdS$d),log=TRUE))
    filterMeans[,k] <- rowMeans(X)  # filter mean
    forecast[,k] <- ym
  }
  new("enkfd_spatPomp", object,
      paramMatrix=array(data=numeric(0),dim=c(0,0)),
      Np=Np,
      filter.mean=filterMeans,
      pred.mean=predMeans,
      forecast=forecast,
      cond.logLik=condlogLik,
      loglik=sum(condlogLik),
      runit_measure = object@runit_measure,
      dunit_measure = object@dunit_measure,
      eunit_measure = object@eunit_measure,
      vunit_measure = object@vunit_measure,
      munit_measure = object@munit_measure,
      unit_names=object@unit_names,
      unit_statenames=object@unit_statenames,
      unit_obsnames = object@unit_obsnames)
}



