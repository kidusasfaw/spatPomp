setGeneric("dunit_measure", function(object,...)standardGeneric("dunit_measure"))

setMethod(
  "dunit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, y, x, unit, time, params, log = TRUE, .gnsi = TRUE, ...){
    pompLoad(object)
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    out <- .Call(do_dunit_measure,object,y,x,time,unit,params,log,.gnsi)
    pompUnload(object)
    out
  }
)
