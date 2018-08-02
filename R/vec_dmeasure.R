## evaluate the unit measurement model density function for each unit

vec_dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  nunits <- length(object@units)
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(params) <- "double"
  retvec <- vector(length = nunits)
  for(i in 1:nunits){
    retvec[i] <- .Call(do_unit_dmeasure,object,y,x,times,i,params,log,statenames,.getnativesymbolinfo)
  }
  pompUnload(object)
  retvec
}

setMethod(
  "vec_dmeasure",
  signature=signature(object="spatpomp"),
  definition=function (object, y, x, times, params, log = FALSE, ...)
    vec_dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
)
