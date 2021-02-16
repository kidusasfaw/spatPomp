#' runit_measure
#'
#' \code{runit_measure} simulates a unit's observation given the entire state
#' @name runit_measure
#' @rdname runit_measure
#' @include spatPomp_class.R spatPomp.R
#' @param object An object of class \code{spatPomp}
#' @param x A state vector for all units
#' @param unit The unit for which to simulate an observation
#' @param time The time for which to simulate an observation
#' @param params parameters to use to simulate an observation
#' @return A matrix with the simulation
#' @examples
#' b <- bm(U=3)
#' s <- states(b)[,1,drop=FALSE]
#' rownames(s) -> rn
#' dim(s) <- c(3,1,1)
#' dimnames(s) <- list(variable=rn, rep=NULL)
#' p <- coef(b); names(p) -> rnp
#' dim(p) <- c(length(p),1); dimnames(p) <- list(param=rnp)
#' o <- obs(b)[,1,drop=FALSE]
#' runit_measure(b, x=s, unit=2, time=1, params=p)
NULL

setGeneric("runit_measure", function(object,...)standardGeneric("runit_measure"))

##' @name runit_measure-spatPomp
##' @aliases runit_measure,spatPomp-method
##' @rdname runit_measure
##' @export
setMethod(
  "runit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, unit, time, params, log = FALSE, .gnsi=TRUE){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    storage.mode(unit) <- "integer"
    out<-.Call(do_runit_measure,
               object,
               x,
               time,
               unit,
               params,
               .gnsi)[unit,,1,drop=FALSE]
    pompUnload(object)
    out
  }
)
