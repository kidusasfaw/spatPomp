##' Show methods
##'
##' Display the object, according to its class.

##' @name show-spatPomp
##' @aliases show,spatPomp-method
##' @rdname show
##' @param object a \code{spatPomp} object
##' @export
setMethod(
  "show",
  signature=signature(object="spatPomp"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)
