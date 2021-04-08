#' eunit_measure
#'
#' \code{eunit_measure} evaluates the expectation of a unit's observation given the entire state
#' @name eunit_measure
#' @rdname eunit_measure
#' @include spatPomp_class.R spatPomp.R
#' @param object An object of class \code{spatPomp}
#' @param x A state vector for all units
#' @param unit The unit for which to evaluate the expectation
#' @param time The time for which to evaluate the expectation
#' @param log logical; should the density be returned on log scale?
#' @param Np numeric; defaults to 1 and the user need not change this
#' @param params parameters at which to evaluate the unit expectation
#' @return A class \sQuote{matrix} with the unit expected observation for spatial unit \code{unit}
#' corresponding to the corresponding states in \code{x}.

#' @examples
#' b <- bm(U=3)
#' s <- states(b)[,1,drop=FALSE]
#' rownames(s) -> rn
#' dim(s) <- c(3,1,1)
#' dimnames(s) <- list(variable=rn, rep=NULL)
#' p <- coef(b); names(p) -> rnp
#' dim(p) <- c(length(p),1); dimnames(p) <- list(param=rnp)
#' o <- obs(b)[,1,drop=FALSE]
#' eunit_measure(b, x=s, unit=2, time=1, params=p)
#'
NULL

setGeneric("eunit_measure", function(object,...)standardGeneric("eunit_measure"))

##' @name eunit_measure-spatPomp
##' @aliases eunit_measure,spatPomp-method
##' @rdname eunit_measure
##' @export
setMethod(
  "eunit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, unit, time, params, Np=1, log=FALSE){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    out <- .Call('do_theta_to_e',
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
