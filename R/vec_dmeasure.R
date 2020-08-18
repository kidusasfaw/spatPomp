##' Evaluate the unit measurement model density function for each unit
##' @include spatPomp_class.R
##' @rdname vec_dmeasure
##'
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
    weights[i,,] <- .Call(do_unit_dmeasure,object,y,x,times,units[i],params,log,.gnsi)
  }
  pompUnload(object)
  return(weights)
}

##' @export
setMethod(
  "vec_dmeasure",
  signature=signature(object="spatPomp"),
  definition=function (object, y, x, units, times, params, log = FALSE, ...){
    if(missing(units)) units <- seq(unit_names(object))
    vec_dmeasure.internal(object=object,y=y,x=x,units=units,times=times,params=params,log=log,...)
  }
)
