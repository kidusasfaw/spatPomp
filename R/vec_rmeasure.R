##' Vector of simulated measurements
##'
##' Simulate from the unit measurement model density function for each unit
##' @param object a \code{spatPomp} object
##' @param x numeric; state at which measurements are simulated
##' @param times numeric; time at which measurements are simulated
##' @param params numeric; parameter set at which measurements are simulated
##' @param \dots additional parameters will be ignored
##' @name vec_rmeasure
##' @include spatPomp_class.R
##' @rdname vec_rmeasure
NULL

setGeneric("vec_rmeasure", function(object,...)standardGeneric("vec_rmeasure"))

##' @name vec_rmeasure-spatPomp
##' @rdname vec_rmeasure
##' @aliases vec_rmeasure,spatPomp-method
##' @export
setMethod(
  "vec_rmeasure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, times, params, ...){
    vec_rmeasure.internal(object=object,x=x,times=times,params=params,...)
  }
)

vec_rmeasure.internal <- function (object, x, times, params, .gnsi = TRUE, ...) {
  pompLoad(object)
  nunits <- length(unit_names(object))
  nparticles <- ncol(x)
  ntimes <- length(times)
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  weights <- array(dim=c(nunits,nparticles,ntimes))

  for(i in 1:nunits){
    # for girf params is not nparams by nparticles. instead it's npars by nunits by nparticles
    if(length(dim(params)) > 2){
      params <- params[,i,]
    }
    weights[i,,] <- .Call(do_runit_measure,object,x,times,i,params,.gnsi)
  }
  pompUnload(object)
  return(weights)
}
