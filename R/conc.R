##' Concatenate spatPomp objects into a listie
##'
##' @description Internal methods to concatenate objects into useful listie.
##' @details Not exported.
##' @name conc
##' @rdname conc
##' @include listie.R
##' @importFrom stats setNames
##' @keywords internal
##'
NULL

setGeneric(
  "conc",
  function (...)
    standardGeneric("conc")
)

setMethod(
  "conc",
  signature=signature(...="ANY"),
  definition=function (...) {
    if (...length()==0L) {
      cls <- "missing"
    } else {
      cls <- unique(vapply(list(...),class,character(1L)))
    }
    pStop_(
      sQuote("c"),
      " is not defined for objects of ",
      ngettext(length(cls),"class ","classes "),
      paste(sQuote(cls),collapse=", "),"."
    )
  }
)

flatten_list <- function (...) {
  unlist(
    lapply(
      list(...),
      \(z) {
        if (is(z,"list")) {
          setNames(as(z,"list"),names(z))
        } else {
          z
        }
      }
    )
  )
}

##' @rdname conc
setMethod(
  "conc",
  signature=signature(...="SpatPomp"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "spatPompList",
        lapply(y,as,"spatPomp")
      ),
      names(y)
    )
  }
)

##' @rdname conc
setMethod(
  "conc",
  signature=signature(...="Bpfilter"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "bpfilterList",
        lapply(y,as,"bpfilterd_spatPomp")
      ),
      names(y)
    )
  }
)

##' @rdname conc
setMethod(
  "conc",
  signature=signature(...="Ibpf"),
  definition=function (...) {
    y <- flatten_list(...)
    setNames(
      new(
        "ibpfList",
        lapply(y,as,"ibpfd_spatPomp")
      ),
      names(y)
    )
  }
)
