##' Coerce to spatPomp
##'
##' Convert to class spatPomp object
##'
##' @name as_spatPomp
##' @rdname as_spatPomp
##' @include spatPomp_class.R
##' @return a class \sQuote{spatPomp} representation of the object.
NULL

##' @name coerce-pomp-spatPomp
##' @aliases coerce,pomp,spatPomp-method
##' @rdname as_spatPomp
##'
##' @details
##' When \code{object} is a simple \sQuote{pomp} object,
##' construct and return a one-dimensional \sQuote{spatPomp} object.
##'
setAs(
  from="pomp",
  to="spatPomp",
  def = function (from) {
    unit_sn <- character(0)
    if(dim(from@states)[1] > 0) unit_sn <- rownames(from@states)
    unit_on <- character(0)
    if(dim(from@data)[1] > 0) unit_on <- rownames(from@data)
    unit_cn <- character(0)
    if(length(get_covariate_names(from@covar))>0) unit_cn <- get_covariate_names(from@covar)

    new("spatPomp",
        from,
        unit_names="unit1",
        unitname="unit",
        unit_statenames=unit_sn,
        unit_obsnames = unit_on,
        unit_covarnames = unit_cn)
  }
)
