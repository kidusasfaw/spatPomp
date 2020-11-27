setGeneric("munit_measure", function(object,...)standardGeneric("munit_measure"))

setMethod(
  "munit_measure",
  signature=signature(object="spatPomp"),
  definition=function (object, x, vc, time, params, Np=1, gnsi=TRUE){
    pompLoad(object)
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    storage.mode(vc) <- "double"
    out <- .Call('do_v_to_theta',
          object=object,
          X=x,
          vc=vc,
          Np = Np,
          times=time,
          params=params,
          gnsi=gnsi)
    pompUnload(object)
    out
  }
)
