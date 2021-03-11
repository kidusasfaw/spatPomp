##' Undefined
##'
##' Check for undefined methods.
##'
##' @return
##' Returns \code{TRUE} if the \pkg{pomp} workhorse method is undefined,
##' \code{FALSE} if it is defined,
##' and \code{NA} if the question is inapplicable.
##'
##' @name undefined
##' @rdname undefined
##' @include pstop.R
##' @keywords internal
##'
##' @param object  object to test.
##'
##' @param ... currently ignored.
##'
NULL
pompfunmode <- list(undef=0L,Rfun=1L,native=2L,regNative=3L)
skeletontype <- list(undef=0L,vectorfield=1L,map=2L)
##' @rdname undefined
setGeneric(
  "undefined",
  function (object, ...)
    standardGeneric("undefined")
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="NULL"),
  definition=function (object, ...) TRUE
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="ANY"),
  definition=function (object, ...) NA
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="missing"),
  definition=function (...) NA
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="pomp_fun"),
  definition=function (object, ...) {
    object@mode == pompfunmode$undef
  }
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="partransPlugin"),
  definition=function (object, ...) {
    undefined(object@to) || undefined(object@from)
  }
)

##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="rprocPlugin"),
  definition=function (object, ...) {
    undefined(object@step.fn) && undefined(object@rate.fn)
  }
)
##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="covartable"),
  definition=function (object) {
    nrow(object@table) == 0L
  }
)

##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="skelPlugin"),
  definition=function (object, ...) {
    object@type==skeletontype$undef || undefined(object@skel.fn)
  }
)
