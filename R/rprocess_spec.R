## also defined in 'pomp_internal.h'
rprocmode <- list(default=0L,onestep=1L,discrete=2L,euler=3L,gillespie=4L)

setClass(
  "rprocPlugin",
  slots=c(
    csnippet='logical',
    slotname='character',
    type='integer',
    step.fn="ANY",
    rate.fn="ANY"
  ),
  prototype=prototype(
    csnippet=FALSE,
    slotname=character(0),
    type=rprocmode$default,
    step.fn=NULL,
    rate.fn=NULL
  )
)

setClass(
  "onestepRprocPlugin",
  contains="rprocPlugin",
  prototype=prototype(
    type=rprocmode$onestep
  )
)

setClass(
  "discreteRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    delta.t="numeric"
  ),
  prototype=prototype(
    type=rprocmode$discrete,
    delta.t=1.0
  )
)

setClass(
  "eulerRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    delta.t="numeric"
  ),
  prototype=prototype(
    type=rprocmode$euler,
    delta.t=NA_real_
  )
)

setClass(
  "gillespieRprocPlugin",
  contains="rprocPlugin",
  slots=c(
    hmax="numeric",
    v="matrix"
  ),
  prototype=prototype(
    type=rprocmode$gillespie
  )
)

rproc_plugin <- function (object, step.fn, rate.fn) {
  if (missing(object)) {
    new("rprocPlugin")
  } else {
    if (!missing(step.fn)) object@step.fn <- step.fn
    if (!missing(rate.fn)) object@rate.fn <- rate.fn
    object
  }
}
