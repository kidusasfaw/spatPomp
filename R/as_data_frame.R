##' Coerce to data frame
##'
##' \pkg{spatPomp} objects can be recast as data frames.
##'
##' @name as.data.frame
##' @rdname as_data_frame
##' @include spatPomp_class.R
##' @aliases coerce,spatPomp,data.frame-method
##' @importFrom rlang `:=` .data
##'
##' @details
##' When \code{object} is a simple \sQuote{spatPomp} object,
##' \code{as(object,"data.frame")} or \code{as.data.frame(object)} results in a
##' data frame with the times, units, observables, states (if known), and
##' interpolated covariates (if any).
##' @param x a \code{spatPomp} object.
setAs(
  from="spatPomp",
  to="data.frame",
  def = function (from) {
    # get all names
    cnames <- pomp:::get_covariate_names(from@covar)
    unitname <- from@unitname
    timename <- from@timename

    # set up unit names for obs, states and covars
    unit_stateobscovars <- c(from@unit_obsnames)

    # get the observation, covars (if any) and states (if any)
    dat <- t(obs(from))
    if (length(from@states)>0) {
      nm <- colnames(dat)
      dat <- cbind(dat,t(from@states))
      colnames(dat) <- c(nm,rownames(from@states))
      unit_stateobscovars <- c(unit_stateobscovars, from@unit_statenames)
    }
    if (length(cnames) > 0) {
      nm <- c(colnames(dat),cnames)
      y <- .Call('lookup_in_table',from@covar,from@times,PACKAGE = 'pomp')
      dat <- cbind(dat,t(y))
      colnames(dat) <- nm
      unit_stateobscovars <- c(unit_stateobscovars, from@unit_covarnames)
    }
    # function to split unit name and unit index
    unit_stateobscovars_pat <- paste0(paste("^",unit_stateobscovars,sep=""), collapse = "|")
    get_unit_index_from_name <- function(name){
      stringr::str_split(name,unit_stateobscovars_pat)[[1]][2]
    }
    get_unit_index_from_name_v <- Vectorize(get_unit_index_from_name)

    # turn into data.frame (from matrix) and complete with time name
    dat <- as.data.frame(dat)
    dat <- cbind(from@times,dat)
    colnames(dat)[1] <- timename
    # convert to long format with column for stateobscovars
    no_time_colnames <- colnames(dat)[-1]
    shared_covnames_ix <- which(no_time_colnames %in% from@shared_covarnames)
    if(length(shared_covnames_ix) > 0)
      to_gather <- no_time_colnames[-shared_covnames_ix]
    else
      to_gather <- no_time_colnames
    to_select <- c(timename, unitname, "stateobscovars", "val")
    to_arrange <- rlang::syms(c(timename, unitname, "stateobscovars"))
    to_final_select <- c(timename, unitname, unit_stateobscovars)
    gathered <- dat %>%
      tidyr::gather_(key="stateobscovars", val="val", to_gather) %>%
      dplyr::mutate(ui = get_unit_index_from_name_v(.data$stateobscovars)) %>%
      dplyr::mutate(!!unitname := unit_names(from)[as.integer(.data$ui)]) %>%
      dplyr::select(-.data$ui) %>%
      dplyr::arrange(!!!to_arrange)

    # get the type of stateobscovars from the stateobscovars column
    stateobscovarstype <- sapply(gathered$stateobscovars,
                           FUN=function(x) stringr::str_extract(
                             x,unit_stateobscovars_pat))
    gathered$stateobscovarstype <- stateobscovarstype

    # spread stateobscovartype column to get columns for all unitnames
    gathered <- gathered %>%
      dplyr::select(-.data$stateobscovars) %>%
      tidyr::spread(key = .data$stateobscovarstype, value = .data$val)%>%
      dplyr::select(to_final_select) %>%
      dplyr::arrange(!!rlang::sym(timename),
                     match(!!rlang::sym(unitname), unit_names(from)))
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

