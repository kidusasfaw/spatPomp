setGeneric("dunit_measure", function(object,...)standardGeneric("dunit_measure"))
#' dunit_measure
#'
#' \code{dunit_measure} evaluates the unit measurement density of a unit's observation given the entire state
#'
#' @name dunit_measure
#' @rdname dunit_measure
#' @include spatPomp_class.R spatPomp.R
#' @param object An object of class \code{spatPomp}
#' @param y A U by 1 matrix of observations for all units
#' @param x A state vector for all units
#' @param unit The unit for which to evaluate the unit measurement density
#' @param params parameters at which to evaluate the unit measurement density
#' @return A matrix with the unit measurement density
#' @examples
#' b <- bm(U=3)
#' s <- states(b)[,1,drop=FALSE]
#' rownames(s) -> rn
#' dim(s) <- c(3,1,1)
#' dimnames(s) <- list(variable=rn, rep=NULL)
#' p <- coef(b); names(p) -> rnp
#' dim(p) <- c(length(p),1); dimnames(p) <- list(param=rnp)
#' o <- obs(b)[,1,drop=FALSE]
#' dunit_measure(b, y=o, x=s, unit=1, time=1, params=p)
#' @export
setMethod(
  "dunit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, y, x, unit, time, params, log = TRUE, .gnsi = TRUE, ...){
    pompLoad(object)
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    storage.mode(unit) <- "integer"
    storage.mode(params) <- "double"
    out <- .Call(do_dunit_measure,object,y,x,time,unit,params,log,.gnsi)
    pompUnload(object)
    out
  }
)
