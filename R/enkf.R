##' Ensemble Kalman filters
##'
##' The ensemble Kalman filter
##'
##' @name senkf
##' @rdname senkf
##' @include spatpomp_class.R spatpomp.R
##' @aliases senkf  senkf,ANY-method senkf,missing-method
##' @author Kidus Asfaw
##' @family particle filtering methods
##' @family \pkg{spatpomp} parameter estimation methods
##'
##' @inheritParams spatpomp
##' @param Np the number of particles to use.
##' @param h function returning the expected value of the observation given the
##' state.
##' @param C matrix converting state vector into expected value of the
##' observation.
##' @param R matrix; variance of the measurement noise.
##'
##' @return
##' An object of class \sQuote{kalmand_spatpomp}.
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
  "kalmand_spatpomp",
  contains="kalmand_pomp",
  slots=c(
    units = 'character',
    unit_index = 'character',
    unit_statenames = 'character',
    global_statenames = 'character',
    obstypes = 'character',
    unit_dmeasure = 'pomp_fun',
    unit_rmeasure = 'pomp_fun'
  ),
  prototype=prototype(
    unit_dmeasure = pomp2:::pomp_fun(slotname="unit_dmeasure"),
    unit_rmeasure = pomp2:::pomp_fun(slotname="unit_rmeasure")
  )
)

setGeneric(
  "senkf",
  function (data, ...)
    standardGeneric("senkf")
)


setMethod(
  "senkf",
  signature=signature(data="missing"),
  definition=function (...) {
    pomp2:::reqd_arg("senkf","data")
  }
)

setMethod(
  "senkf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("senkf",data)
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

##' @name senkf-spatpomp
##' @aliases senkf,spatpomp-method
##' @rdname senkf
##' @export
setMethod(
  "senkf",
  signature=signature(data="spatpomp"),
  function (data,
    Np, h, R,
    ..., verbose = getOption("verbose", FALSE)) {
    tryCatch(
      kp <- pomp2:::enkf.internal(
        data,
        Np=Np,
        h=h,
        R=R,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp2:::pStop("senkf",conditionMessage(e))
    )
    new("kalmand_spatpomp", kp,
        unit_rmeasure = data@unit_rmeasure,
        unit_dmeasure = data@unit_dmeasure,
        units=data@units,
        unit_index=data@unit_index,
        unit_statenames=data@unit_statenames,
        global_statenames=data@global_statenames,
        obstypes = data@obstypes)

  }
)
