#' @include spatpomp_class.R
#'
## this file contains some basic methods definitions
setAs(
  from="pomp",
  to="spatpomp",
  def = function (from) {
    new("spatpomp",from,
        unit_dmeasure=from@dmeasure,
        units="unit",
        unit_index="unit",
        unit_statenames=character(0),
        global_statenames=character(0),
        obstypes = rownames(from@data))
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
  "print",
  signature=signature(x="spatpomp"),
  definition=function (x, ...) {
    cat("<object of class ",sQuote("spatpomp"),">\n",sep="")
    invisible(x)
  }
)

setMethod(
  "simulate",
  signature=signature(object="spatpomp"),
  definition=function(object, nsim = 1, seed = NULL,
                       format = c("spatpomps", "arrays", "data.frame"),
                       include.data = FALSE,...) {
    s1 <- simulate(pomp(object), ...)
    if(format=="data.frame"){
      # convert to long format and output
    }
    if(format=="spatpomps"){
      # add back spatpomp components
    }
    s1
  }
)

setMethod(
  "show",
  signature=signature(object="spatpomp"),
  definition=function (object) {
    cat(length(object@times),"records of",
        nrow(obs(object)),"observables,",
        "recorded from t =",
        min(object@times),"to",max(object@times), "for ", length(object@units), " units\n")
    cat("summary of data:\n")
    print(summary(as.data.frame(t(obs(object)))))
    cat("zero time, t0 = ",object@t0,"\n",sep="")
    if (!is.null(dim(object@covar))) {
      cat(nrow(object@covar),"records of",
          ncol(object@covar),"covariates,",
          "recorded from t =",min(object@tcovar),
          "to",max(object@tcovar),"\n")
      cat("summary of covariates:\n")
      print(summary(as.data.frame(object@covar)))
    }
    cat("process model simulator, rprocess = ")
    show(object@rprocess)
    cat("process model density, dprocess = ")
    show(object@dprocess)
    cat("measurement model simulator, rmeasure = ")
    show(object@rmeasure)
    cat("measurement model density, dmeasure = ")
    show(object@dmeasure)
    cat("prior simulator, rprior = ")
    show(object@rprior)
    cat("prior density, dprior = ")
    show(object@dprior)
    show(object@skeleton)
    cat("rinit = ")
    show(object@rinit)
    cat("parameter transformation  = ")
    show(object@partrans)
    if (length(coef(object))>0) {
      cat("parameter(s):\n")
      print(coef(object))
    } else {
      cat ("parameter(s) unspecified\n");
    }
    if (length(object@userdata)>0) {
      cat("extra user-defined variables: ",
          paste(sapply(names(object@userdata),sQuote),collapse=", "),
          "\n")
    }
    invisible(NULL)
  }
)

