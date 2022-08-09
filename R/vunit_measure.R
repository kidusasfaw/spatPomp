#' vunit_measure
#'
#' \code{vunit_measure} evaluates the variance of a unit's observation given the entire state
#' @name vunit_measure
#' @rdname vunit_measure
#' @include spatPomp_class.R spatPomp.R
#' @param object An object of class \code{spatPomp}
#' @param x A state vector for all units
#' @param unit The unit for which to evaluate the variance
#' @param time The time for which to evaluate the variance
#' @param Np numeric; defaults to 1 and the user need not change this
#' @param params parameters at which to evaluate the unit variance
#' @return A matrix with the unit measurement variance implied by the state, \code{x},
#' and the parameter set \code{params} for unit \code{unit}.
#' @examples Complete examples are provided in the package tests
#' \dontrun{
#' b <- bm(U=3)
#' s <- states(b)[,1,drop=FALSE]
#' rownames(s) -> rn
#' dim(s) <- c(3,1,1)
#' dimnames(s) <- list(variable=rn, rep=NULL)
#' p <- coef(b); names(p) -> rnp
#' dim(p) <- c(length(p),1); dimnames(p) <- list(param=rnp)
#' o <- obs(b)[,1,drop=FALSE]
#' vunit_measure(b, x=s, unit=2, time=1, params=p)
#' }
NULL
setGeneric("vunit_measure", function(object,...)standardGeneric("vunit_measure"))

##' @name vunit_measure-spatPomp
##' @aliases vunit_measure,spatPomp-method
##' @rdname vunit_measure
##' @export
setMethod(
  "vunit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, unit, time, params, Np=1){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    out <- .Call(do_theta_to_v,
          object=object,
          X=x,
          Np = as.integer(Np),
          times=time,
          params=params,
          gnsi=TRUE)
    pompUnload(object)
    out[unit,,,drop=FALSE]
  }
)
