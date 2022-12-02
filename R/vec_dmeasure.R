##' Vector of measurement densities
##'
##' Evaluate the unit measurement model density function for each unit.
##' This method is used primarily as part of likelihood evaluation and parameter inference algorithms.
##' @param object a \code{spatPomp} object
##' @param y numeric; measurements whose densities given the latent states are evaluated
##' @param x numeric; state at which conditional measurement densities are evaluated
##' @param units numeric; units at which measurement densities are evaluated
##' @param times numeric; time at which measurement densities are evaluated
##' @param params numeric; parameter set at which measurement densities is evaluated
##' @param log logical; should the outputted measurement densities be on log scale?
##' @param \dots additional parameters will be ignored
##' @name vec_dmeasure
##' @include spatPomp_class.R
##' @rdname vec_dmeasure
##' @author Kidus Asfaw
##' @return An array of dimension \code{length(unit_names(object))} by \code{dim(x)[2]} by \code{dim(x)[3]}
##' representing each unit's measurement density assessed for each replicate in \code{x} for each observation time.
NULL

setGeneric("vec_dmeasure", function(object,...)standardGeneric("vec_dmeasure"))

##' @name vec_dmeasure-spatPomp
##' @rdname vec_dmeasure
##' @aliases vec_dmeasure,spatPomp-method
##' @export
setMethod(
  "vec_dmeasure",
  signature=signature(object="spatPomp"),
  definition=function (object, y, x, units, times, params, log = FALSE, ...){
    if(missing(units)) units <- seq(unit_names(object))
    vec_dmeasure.internal(object=object,y=y,x=x,units=units,times=times,params=params,log=log,...)
  }
)

vec_dmeasure.internal <- function (object, y, x, units, times, params, log = FALSE, .gnsi = TRUE, ...) {
  pompLoad(object)
  nunits <- length(units)
  nparticles <- ncol(x)
  ntimes <- length(times)
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  weights <- array(dim=c(nunits,nparticles,ntimes))

  for(i in seq(nunits)){
    # for girf params is not nparams by nparticles. instead it's npars by nunits by nparticles
    if(length(dim(params)) > 2){
      params <- params[,i,]
    }
    weights[i,,] <- .Call(do_dunit_measure,object,y,x,times,
      as.integer(units[i]-1),params,log,.gnsi)
  }
  pompUnload(object)
  return(weights)
}
