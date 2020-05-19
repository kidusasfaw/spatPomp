##' Ensemble Kalman filters
##'
##' A function to perform filtering using the ensemble Kalman filter of Evensen, G. (1994)
##'
##' @name senkf2
##' @rdname senkf2
##' @include spatPomp_class.R spatPomp.R
##' @aliases senkf2  senkf2,ANY-method senkf2,missing-method
##' @family particle filtering methods
##' @family \pkg{spatPomp} parameter estimation methods
##'
##' @inheritParams spatPomp
##' @param Np the number of particles to use.
##'
##' @return
##' An object of class \sQuote{senkfd_spatPomp}.
##'
##' @references
##' Evensen, G. (1994) Sequential data assimilation with a
##' nonlinear quasi-geostrophic model using Monte Carlo methods to forecast
##' error statistics Journal of Geophysical Research: Oceans 99:10143--10162
##'
##' Evensen, G. (2009) Data assimilation: the ensemble Kalman filter
##' Springer-Verlag.
##'
##' Anderson, J. L. (2001) An Ensemble Adjustment Kalman Filter for Data
##' Assimilation Monthly Weather Review 129:2884--2903
NULL

setClass(
  "senkfd_spatPomp",
  contains="kalmand_pomp",
  slots=c(
    units = 'character',
    unit_statenames = 'character',
    obstypes = 'character',
    unit_dmeasure = 'pomp_fun',
    unit_rmeasure = 'pomp_fun',
    unit_emeasure = 'pomp_fun',
    unit_vmeasure = 'pomp_fun',
    unit_mmeasure = 'pomp_fun'
  ),
  prototype=prototype(
    unit_dmeasure = pomp:::pomp_fun(slotname="unit_dmeasure"),
    unit_rmeasure = pomp:::pomp_fun(slotname="unit_rmeasure"),
    unit_emeasure = pomp:::pomp_fun(slotname="unit_emeasure"),
    unit_vmeasure = pomp:::pomp_fun(slotname="unit_vmeasure"),
    unit_mmeasure = pomp:::pomp_fun(slotname="unit_mmeasure")
  )
)

setMethod(
  "senkf2",
  signature=signature(data="missing"),
  definition=function (...) {
    pomp:::reqd_arg("senkf2","data")
  }
)

setMethod(
  "senkf2",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("senkf2",data)
  }
)

## ENSEMBLE KALMAN FILTER (ENKF)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$
## Forecast: $Y_t=h(X_t)$
## Forecast mean: $N_t=\langle Y \rangle$.
## Forecast variance: $S_t=\langle\langle Y \rangle\rangle$
## State/forecast covariance: $W_t=\langle\langle X,Y\rangle\rangle$
## Kalman gain: $K_t = W_t\,S_t^{-1}$
## New observation: $y_t\in \mathbb{R}^{n\times 1}$
## Updated ensemble: $X^u_{t}=X_t + K_t\,(O_t - Y_t)$
## Filter mean: $m_t=\langle X^u_t \rangle = \frac{1}{q} \sum\limits_{i=1}^q x^{u_i}_t$

##' @name senkf-spatPomp
##' @aliases senkf,spatPomp-method
##' @rdname senkf
##' @export
setMethod(
  "senkf2",
  signature=signature(data="spatPomp"),
  function (data,
            Np,
            ..., verbose = getOption("verbose", FALSE)) {
    tryCatch(
      sp <- senkf.internal(
        data,
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp:::pStop("senkf2",conditionMessage(e))
    )
  }
)

senkf.internal <- function (object,
                           Np,
                           ..., verbose) {

  verbose <- as.logical(verbose)

  # object <- pomp(object,...,verbose=verbose)

  if (pomp:::undefined(object@rprocess))
    pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")
  if (pomp:::undefined(object@unit_emeasure))
    pStop_(paste(sQuote(c("unit_emeasure")),collapse=", ")," is a needed basic component.")
  if (pomp:::undefined(object@unit_vmeasure))
    pStop_(paste(sQuote(c("unit_vmeasure")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  #R <- as.matrix(R)
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

    # ensemble of data
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
    meas_var <- tryCatch(
      .Call('do_theta_to_v',
            object=object,
            X=X,
            Np = as.integer(Np[1]),
            times=tt[k+1],
            params=params,
            gnsi=TRUE),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE) # nocov
      }
    )
    dim(meas_var) <- c(length(spat_units(object)),  Np[1])
    R <- diag(rowMeans(meas_var))
    sqrtR <- tryCatch(
      t(chol(R)),                     # t(sqrtR)%*%sqrtR == R
      error = function (e) {
        pStop_("degenerate ",sQuote("R"), "at time ", squote(k), ": ",conditionMessage(e))
      }
    )
    X <- X[,,1]
    predMeans[,k] <- pm <- rowMeans(X) # prediction mean

    dim(Y) <- c(length(spat_units(object)), Np[1])


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

    condlogLik[k] <- sum(dnorm(x=crossprod(svdS$u,resid2),mean=0,sd=sqrt(svdS$d),log=TRUE))
    filterMeans[,k] <- rowMeans(X)  # filter mean
    forecast[,k] <- ym
  }
  new("senkfd_spatPomp", object, Np=Np,
      filter.mean=filterMeans,
      pred.mean=predMeans,
      forecast=forecast,
      cond.loglik=condlogLik,
      loglik=sum(condlogLik),
      unit_rmeasure = object@unit_rmeasure,
      unit_dmeasure = object@unit_dmeasure,
      unit_emeasure = object@unit_emeasure,
      unit_vmeasure = object@unit_vmeasure,
      unit_mmeasure = object@unit_mmeasure,
      units=object@units,
      unit_statenames=object@unit_statenames,
      obstypes = object@obstypes)
}



