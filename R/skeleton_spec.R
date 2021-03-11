skeletontype <- list(undef=0L,vectorfield=1L,map=2L)

setClass(
  "skelPlugin",
  slots=c(
    csnippet='logical',
    type='integer',
    skel.fn="ANY"
  ),
  prototype=prototype(
    csnippet=FALSE,
    type=skeletontype$undef,
    skel.fn=NULL
  )
)

setClass(
  "vectorfieldPlugin",
  contains="skelPlugin",
  prototype=prototype(
    type=skeletontype$vectorfield
  )
)

setClass(
  "mapPlugin",
  contains="skelPlugin",
  slots=c(
    delta.t="numeric"
  ),
  prototype=prototype(
    type=skeletontype$map,
    delta.t=1.0
  )
)

skel_plugin <- function (object, skel.fn) {
  if (missing(object)) {
    new("skelPlugin")
  } else {
    if (!missing(skel.fn)) object@skel.fn <- skel.fn
    object
  }
}
