setGeneric("runit_measure", function(object,...)standardGeneric("runit_measure"))

setMethod(
  "runit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, unit, time, params, log = FALSE, .gnsi=TRUE){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    out<-.Call(do_runit_measure,
               object,
               x,
               time,
               unit,
               params,
               .gnsi)
    pompUnload(object)
    out
  }
)
