##' Unit names of a spatiotemporal model
##'
##' \code{unit_names} outputs the contents of the \code{unit_names} slot
##' of a \code{spatPomp} object. The order in which the units
##' appear in the output vector determines the order in which latent
##' states and observations for the spatial units are stored.
##'
##' @name unit_names
##' @rdname unit_names
##' @include spatPomp_class.R
##' @return A character vector with the unit names used to create the \sQuote{spatPomp} object.
NULL
setGeneric("unit_names", function(x)standardGeneric("unit_names"))

##' @name unit_names-spatPomp
##' @rdname unit_names
##' @aliases unit_names,spatPomp-method
##' @param x a \code{spatPomp} object
##' @export
setMethod(
  "unit_names",
  signature=signature(x="spatPomp"),
  definition=function(x) x@unit_names
)
