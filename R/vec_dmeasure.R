## evaluate the unit measurement model density function for each unit

vec_dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  nunits <- length(object@units)
  nparticles <- ncol(x)
  ntimes <- length(times)
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  weights <- array(dim=c(nunits,nparticles,ntimes))

  for(i in 1:nunits){
    # get relevant observation types and state types
    # for(i in object@obstypes){
    #   relevant_y <- y[paste0(object@obstypes,i),,drop=FALSE]
    #   row.names(relevant_y) <- object@obstypes
    # }
    print('vec_dmeasure is not returning for all units')
    return(.Call(do_unit_dmeasure,object,y,x,times,i,params,log,.getnativesymbolinfo))
    #weights[i,,] <- .Call(do_unit_dmeasure,object,y,x,times,i,params,log,.getnativesymbolinfo)
  }
  pompUnload(object)
  return(weights)
}

setMethod(
  "vec_dmeasure",
  signature=signature(object="spatpomp"),
  definition=function (object, y, x, times, params, log = FALSE, ...)
    vec_dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
)
