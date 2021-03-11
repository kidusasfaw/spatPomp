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
##' @name show-skelPlugin
##' @aliases show,skelPlugin-method
##' @rdname show
##' @param object a \code{skelPlugin} object
##' @export
setMethod(
  "show",
  signature=signature(object="skelPlugin"),
  definition=function (object) {
    cat("<default>\n\n")
  }
)
##' @name show-vectorfieldPlugin
##' @aliases show,vectorfieldPlugin-method
##' @rdname show
##' @param object a \code{vectorfieldPlugin} object
##' @export
setMethod(
  "show",
  signature=signature(object="vectorfieldPlugin"),
  definition=function (object) {
    cat("vectorfield:\n  - ")
    show(object@skel.fn)
  }
)
##' @name show-mapPlugin
##' @aliases show,mapPlugin-method
##' @rdname show
##' @param object a \code{mapPlugin} object
##' @export
setMethod(
  "show",
  signature=signature(object="mapPlugin"),
  definition=function (object) {
    cat("map:\n")
    cat("  - time-step =",object@delta.t,"\n")
    cat("  - ")
    show(object@skel.fn)
  }
)
