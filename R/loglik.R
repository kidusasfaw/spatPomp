##' Log likelihood
##'
##' Extract the estimated log likelihood.
##'
##' @name logLik
##' @rdname loglik
##' @aliases logLik logLik,ANY-method logLik,missing-method
##' @include girf.R bpfilter.R abf.R iubf.R abfir.R igirf.R
##'
##' @param object fitted model object
##'
##' @return a numeric which is an estimated log likelihood
NULL

##' @name logLik-girfd_spatPomp
##' @aliases logLik,girfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="girfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-bpfilterd_spatPomp
##' @aliases logLik,bpfilterd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="bpfilterd_spatPomp"),
  definition=function(object)object@loglik
)


##' @name logLik-abfd_spatPomp
##' @aliases logLik,abfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="abfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-iabf6d_spatPomp
##' @aliases logLik,iabf6d_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="iubfd_spatPomp"),
  definition=function(object)object@loglik
)


##' @name logLik-abfird_spatPomp
##' @aliases logLik,abfird_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="abfird_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-igirfd_spatPomp
##' @aliases logLik,igirfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="igirfd_spatPomp"),
  definition=function(object)object@loglik
)
