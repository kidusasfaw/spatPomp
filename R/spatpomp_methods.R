#' @include spatpomp_class.R
#'
## this file contains some basic methods definitions

## 'coerce' method: allows for coercion of a "pomp" object to a data-frame
setAs(
  from="spatpomp",
  to="data.frame",
  def = function (from) {
    x <- as.data.frame(cbind(from@times,t(from@data)))
    names(x) <- c("time",rownames(from@data))
    if (length(from@states)>0) {
      nm <- names(x)
      x <- cbind(x,t(from@states))
      names(x) <- c(nm,rownames(from@states))
    }
    if (length(from@covar)>0) {
      nm <- c(names(x),colnames(from@covar))
      y <- .Call(lookup_in_table,from@tcovar,from@covar,from@times)
      x <- cbind(x,t(y))
      names(x) <- nm
    }
    x
  }
)

as.data.frame.pomp <- function (x, row.names, optional, ...) as(x,"data.frame")

## parameter transformations
partrans.internal <- function (object, params,
                               dir = c("fromEstimationScale",
                                       "toEstimationScale",
                                       "forward","inverse"),
                               .getnativesymbolinfo = TRUE, ...) {
  if (!object@has.trans) return(params)
  pompLoad(object)
  dir <- switch(match.arg(dir),fromEstimationScale=1L,toEstimationScale=-1L,
                forward=1L,inverse=-1L)
  rv <- .Call(do_partrans,object,params,dir,.getnativesymbolinfo)
  pompUnload(object)
  rv
}

setMethod(
  "partrans",
  signature=signature(object="spatpomp"),
  definition=function (object, params, dir = c("fromEstimationScale",
                                               "toEstimationScale", "forward","inverse"),
                       ...)
    partrans.internal(object=object,params=params,dir=dir,...)
)


obs.internal <- function (object, vars, ...) {
  varnames <- colnames(object@data)
  if (missing(vars))
    vars <- varnames
  else if (!all(vars%in%varnames))
    stop("in ",sQuote("obs"),": some elements of ",
         sQuote("vars")," correspond to no observed variable",call.=FALSE)
  y <- object@data[,vars,,drop=FALSE]
  dimnames(y) <- list(unit = unit(object), variable=colnames(y),time=time(object))
  y
}

## a simple method to extract the data array
setMethod(
  "obs",
  signature=signature(object="spatpomp"),
  definition=obs.internal
)

## a simple method to extract the array of states
states.internal <- function (object, vars, ...) {
  if (length(object@states)==0) {
    NULL
  } else {
    if (missing(vars))
      vars <- seq(length=nrow(object@states))
    x <- object@states[vars,,drop=FALSE]
    dimnames(x) <- list(variable=rownames(x),time=time(object))
    x
  }
}

setMethod(
  "states",
  signature=signature(object="spatpomp"),
  definition=states.internal
)

## a simple method to extract the vector of times
setMethod(
  "time",
  signature=signature(x="spatpomp"),
  definition=function (x, t0 = FALSE, ...) {
    if (t0) c(x@t0,x@times) else x@times
  }
)

## modify the time
setMethod(
  "time<-",
  signature=signature(object="spatpomp"),
  definition=function (object, t0 = FALSE, ..., value) {
    ep <- paste0("in ",sQuote("time<-"),": ")
    if (!is.numeric(value))
      stop(ep,sQuote("value")," must be a numeric vector",call.=FALSE)
    storage.mode(value) <- "double"
    tt <- object@times
    dd <- object@data
    ss <- object@states
    if (t0) {
      object@t0 <- value[1]
      object@times <- value[-1]
    } else {
      object@times <- value
    }
    if (!all(diff(object@times)>0))
      stop(ep,"the times specified must be an increasing sequence",call.=FALSE)
    if (object@t0>object@times[1])
      stop(ep,"the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=FALSE)
    object@data <- array(
      data=NA,
      dim=c(nrow(dd),length(object@times)),
      dimnames=list(rownames(dd),NULL)
    )
    object@data[,object@times%in%tt] <- dd[,tt%in%object@times]
    if (length(ss)>0) {
      object@states <- array(
        data=NA,
        dim=c(nrow(ss),length(object@times)),
        dimnames=list(rownames(ss),NULL)
      )
      for (kt in seq_along(object@times)) {
        wr <- which(object@times[kt]==tt)
        if (length(wr)>0)
          object@states[,kt] <- ss[,wr[1]]
      }
    }
    object
  }
)

## extract the vector of units
setMethod(
  "unit",
  signature=signature(x="spatpomp"),
  definition=function(x,...) x@units
)

## extract the unit index
setMethod(
  "unit_ix",
  signature=signature(x="spatpomp"),
  definition=function(x,...) x@unit_index
)

setMethod(
  "window",
  signature=signature(x="spatpomp"),
  definition=function (x, start, end, ...) {
    tm <- time(x,t0=FALSE)
    if (missing(start))
      start <- tm[1]
    if (missing(end))
      end <- tm[length(tm)]
    tm <- tm[(tm>=start)&(tm<=end)]
    time(x,t0=FALSE) <- tm
    x
  }
)

setMethod(
  "timezero",
  signature=signature(object="spatpomp"),
  definition=function(object,...)object@t0
)

setMethod(
  "timezero<-",
  signature=signature(object="spatpomp"),
  definition=function(object,...,value) {
    ep <- paste0("in ",sQuote("timezero<-"),": ")
    if (value>object@times[1])
      if (!is.numeric(value) || length(value) > 1)
        stop(ep,"the zero-time ",sQuote("t0"),
             " must be a single number",call.=FALSE)
    if (value > object@times[1])
      stop(ep,"the zero-time ",sQuote("t0"),
           " must occur no later than the first observation",call.=FALSE)
    storage.mode(value) <- "double"
    object@t0 <- value
    object
  }
)



setMethod(
  "print",
  signature=signature(x="spatpomp"),
  definition=function (x, ...) {
    cat("<object of class ",sQuote("pomp"),">\n",sep="")
    invisible(x)
  }
)

setMethod(
  "show",
  signature=signature(object="spatpomp"),
  definition=function (object) {print("donezo")}
)
