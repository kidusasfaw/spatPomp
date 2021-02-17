##' Coerce to spatPomp
##'
##' Convert to class spatPomp object
##'
##' @name as_spatPomp
##' @rdname as_spatPomp
##' @include spatPomp_class.R
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
    new("spatPomp",from,
        dunit_measure=from@dmeasure,
        unit_names="unit",
        unit_statenames=character(0),
        unit_obsnames = rownames(from@data))
  }
)
