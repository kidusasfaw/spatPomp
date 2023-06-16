##' listie
##'
##' List-like objects.
##'
##' @name listie
##' @rdname listie
##' @keywords internal
##' @include spatPomp_class.R
##' @include ibpf.R bpfilter.R
NULL

setClass(
  "spatPompList",
  contains="list",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"spatPomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "bpfilterList",
  contains="spatPompList",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"bpfilterd_spatPomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClass(
  "ibpfList",
  contains="bpfilterList",
  validity=function (object) {
    if (length(object) > 0L) {
      if (!all(vapply(object,is,logical(1L),"ibpfd_spatPomp"))) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": dissimilar objects cannot be combined"
        )
        return(retval)
      }
      d <- sapply(object,\(x)dim(x@traces))
      if (!all(apply(d,1L,diff)==0L)) {
        retval <- paste0(
          "error in ",sQuote("c"),
          ": to be combined, ",sQuote("ibpfd_spatPomp"),
          " objects must have chains of equal length"
        )
        return(retval)
      }
    }
    TRUE
  }
)

setClassUnion("SpatPomp",c("spatPomp","spatPompList"))

setClassUnion("Bpfilter",
  c("bpfilterd_spatPomp","bpfilterList")
)

setClassUnion("Ibpf",c("ibpfd_spatPomp","ibpfList"))

setClassUnion("listie",
  members=c("spatPompList","ibpfList","bpfilterList")
)
