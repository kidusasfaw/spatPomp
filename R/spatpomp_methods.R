#' @include spatpomp_class.R
#'
## this file contains some basic methods definitions

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
    if (length(object@tcovar)>0) {
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
    cat("skeleton ",
        if (object@skeleton.type!="undef")
          paste0("(",object@skeleton.type,") ")
        else "",
        "= ",sep="")
    show(object@skeleton)
    cat("initializer = ")
    show(object@initializer)
    cat("parameter transformation (to estimation scale) = ")
    show(object@to.trans)
    cat("parameter transformation (from estimation scale) = ")
    show(object@from.trans)
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

