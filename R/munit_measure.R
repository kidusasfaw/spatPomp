#' munit_measure
#'
#' \code{munit_measure} returns a moment-matched parameter set given an empirically calculated measurement variance and latent states.
#' This is used in \code{girf()} and \code{igirf()} when they are run with \code{kind='moment'}.
#' @name munit_measure
#' @rdname munit_measure
#' @include spatPomp_class.R spatPomp.R
#' @param object An object of class \code{spatPomp}
#' @param x A state vector for all units
#' @param vc The empirically calculated variance used to perform moment-matching
#' @param Np Number of particle replicates for which to get parameter sets
#' @param unit The unit for which to obtain a moment-matched parameter set
#' @param time The time for which to obtain a moment-matched parameter set
#' @param params parameters to use to obtain a moment-matched parameter set
#' @return An array with dimensions \code{dim(array.params)[1]} by \code{dim(x)[2]} by \code{length(unit)} by\code{length(time)}
#' representing the moment-matched parameter set(s) corresponding to the variance of the measurements, \code{vc}, and the states, \code{x}.
#' @examples
#' b <- bm(U=3)
#' s <- states(b)[,1,drop=FALSE]
#' rownames(s) -> rn
#' dim(s) <- c(3,1,1)
#' dimnames(s) <- list(variable=rn, rep=NULL)
#' p <- coef(b); names(p) -> rnp
#' dim(p) <- c(length(p),1); dimnames(p) <- list(param=rnp)
#' o <- obs(b)[,1,drop=FALSE]
#' array.params <- array(p,
#'                       dim = c(length(p),
#'                               length(unit_names(b)), 1, 1),
#'                       dimnames = list(params = rownames(p)))
#' vc <- c(4, 9, 16); dim(vc) <- c(length(vc), 1, 1)
#' munit_measure(b, x=s, vc=vc, Np=1, unit = 1, time=1, params=array.params)
#'
NULL

setGeneric("munit_measure", function(object,...)standardGeneric("munit_measure"))

##' @name munit_measure-spatPomp
##' @aliases munit_measure,spatPomp-method
##' @rdname munit_measure
##' @export
setMethod(
  "munit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, vc, unit, time, params, Np=1){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    storage.mode(vc) <- "double"
    storage.mode(unit) <- "integer"
    out <- .Call(do_v_to_theta,
          object=object,
          X=x,
          vc=vc,
          Np = Np,
          times=time,
          params=params,
          gnsi=TRUE)[,unit,,,drop=FALSE]
    pompUnload(object)
    out
  }
)
