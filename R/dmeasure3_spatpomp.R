#' @include generics.R
#' @include spatpomp_class.R
#'
## evaluate the measurement model density function

dmeasure.internal3 <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    pompLoad(object)
    rv <- .Call(do_dmeasure3,object,y,x,times,params,log,.getnativesymbolinfo)
    pompUnload(object)
    rv
}

#' @exportMethod
setMethod(
    "dmeasure3",
    signature=signature(object="spatpomp"),
    definition=function (object, y, x, times, params, log = FALSE, ...)
        dmeasure.internal3(object=object,y=y,x=x,times=times,params=params,log=log,...)
)
