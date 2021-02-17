##' Print methods
##'
##' Prints its argument.
##'
##' @name print
##' @rdname print
NULL

setGeneric(
  "print",
  function (x)
    standardGeneric("print")
)

##' @name print-spatPomp
##' @aliases print,spatPomp-method
##' @rdname print
##' @param x a \code{spatPomp} object
##' @export
setMethod(
  "print",
  signature=signature(x="spatPomp"),
  definition=function (x) {
    cat("<object of class ",sQuote("spatPomp"),">\n",sep="")
    invisible(x)
  }
)
