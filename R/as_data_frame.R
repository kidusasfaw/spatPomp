##' Coerce to data frame
##'
##' \pkg{spatPomp} objects can be recast as data frames.
##'
##' @name as.data.frame
##' @rdname as_data_frame
##' @include spatPomp_class.R
NULL

##' @name coerce-spatPomp-data.frame
##' @aliases coerce,spatPomp,data.frame-method
##' @rdname as_data_frame
##'
##' @details
##' When \code{object} is a simple \sQuote{spatPomp} object,
##' \code{as(object,"data.frame")} or \code{as.data.frame(object)} results in a
##' data frame with the times, units, observables, states (if known), and
##' interpolated covariates (if any).
##'
setAs(
  from="spatPomp",
  to="data.frame",
  def = function (from) {
    dat_cov <- obs(from)
    cnames <- pomp:::get_covariate_names(from@covar)
    if (length(cnames) > 0) {
      nm <- c(rownames(dat_cov),cnames)  # perhaps not strictly necessary (but see issue #56)
      y <- .Call('lookup_in_table',from@covar,from@times,PACKAGE = 'pomp')
      dat_cov <- cbind(t(dat_cov),t(y))
      colnames(dat_cov) <- nm
    }
    if (length(from@states)>0) {
      nm <- colnames(dat_cov)
      dat_cov <- cbind(dat_cov,t(from@states))
      colnames(dat_cov) <- c(nm,rownames(from@states))
    }
    dat_cov <- as.data.frame(dat_cov)
    dat_cov <- cbind(from@times,dat_cov)
    colnames(dat_cov)[1] <- from@timename

    unit_stateobs <- c(from@obstypes, from@unit_statenames, from@unit_covarnames)
    unit_stateobs_pat <- paste0(paste("^",unit_stateobs,sep=""), collapse = "|")
    get_unit_index_from_statename <- function(statename){
      stringr::str_split(statename,unit_stateobs_pat)[[1]][2]
    }
    get_unit_index_from_statename_v <- Vectorize(get_unit_index_from_statename)

    # convert to long format and output
    to_gather <- colnames(dat_cov)[2:length(colnames(dat_cov))][!c(colnames(dat_cov)[2:length(colnames(dat_cov))]%in%from@shared_covarnames)]
    to_select <- c(from@timename, "unit", "stateobs", "val")
    to_arrange <- rlang::syms(c(from@timename, "unit", "stateobs"))
    to_final_select <- c(from@timename,"unit",unit_stateobs)
    gathered <- dat_cov %>%
      tidyr::gather_(key="stateobs", val="val", to_gather) %>%
      dplyr::mutate(ui = get_unit_index_from_statename_v(stateobs)) %>%
      dplyr::mutate(unit = spat_units(from)[as.integer(ui)]) %>%
      dplyr::select(to_select) %>%
      dplyr::arrange(!!!to_arrange)

    stateobstype <- sapply(gathered$stateobs,FUN=function(x) stringr::str_extract(x,unit_stateobs_pat))
    gathered$stateobstype <- stateobstype

    gathered <- gathered %>%
      dplyr::select(-stateobs) %>%
      tidyr::spread(key = stateobstype, value = val)%>%
      dplyr::select(to_final_select) %>%
      dplyr::arrange(!!rlang::sym(from@timename), match(unit,spat_units(from)))
    gathered
  }
)

##' @method as.data.frame spatPomp
##' @rdname as_data_frame
##'
##' @inheritParams base::as.data.frame
##' @export
##'
as.data.frame.spatPomp <- function (x, ...) as(x,"data.frame")

