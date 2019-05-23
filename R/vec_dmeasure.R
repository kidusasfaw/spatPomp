#' @include spatPomp_class.R
## evaluate the unit measurement model density function for each unit

vec_dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .gnsi = TRUE, ...) {
  pompLoad(object)
  nunits <- length(object@units)
  nparticles <- ncol(x)
  ntimes <- length(times)
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  weights <- array(dim=c(nunits,nparticles,ntimes))

  for(i in 1:nunits){
    # for girf params is not nparams by nparticles. instead it's npars by nunits by nparticles
    if(length(dim(params)) > 2){
      params <- params[,i,]
    }
    weights[i,,] <- .Call(do_unit_dmeasure,object,y,x,times,i,params,log,.gnsi)
  }
  pompUnload(object)
  return(weights)
}

setMethod(
  "vec_dmeasure",
  signature=signature(object="spatPomp"),
  definition=function (object, y, x, times, params, log = FALSE, ...)
    vec_dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
)
