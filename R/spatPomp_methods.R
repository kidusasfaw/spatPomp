##' Basic methods for the spatPomp class
##'
##' @name spatpomp-methods
##' @title spatpomp-methods
##' @rdname spatPomp-methods
##'
##' @include spatPomp_class.R generics.R
##'
NULL

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

##' @export
setMethod(
  "unit_names",
  signature=signature(x="spatPomp"),
  definition=function(x,...) x@unit_names
)


##' @export
setMethod(
  "print",
  signature=signature(x="spatPomp"),
  definition=function (x, ...) {
    cat("<object of class ",sQuote("spatPomp"),">\n",sep="")
    invisible(x)
  }
)

##' @export
setMethod(
  "simulate",
  signature=signature(object="spatPomp"),
  definition=function(object, nsim = 1, seed = NULL,
                       format = c("spatPomps", "arrays", "data.frame"),
                       include.data = FALSE,...) {
    format <- match.arg(format)
    if(format == 'spatPomps') sims <- simulate(pomp(object), format = 'pomps', nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format == 'data.frame') sims <- simulate(pomp(object), format = format, nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format=="data.frame"){
      unitname <- object@unitname
      unit_stateobs <- c(object@unit_obsnames, object@unit_statenames)
      unit_stateobs_pat <- paste0(paste("^",unit_stateobs,sep=""), collapse = "|")
      get_unit_index_from_statename <- function(statename){
        stringr::str_split(statename,unit_stateobs_pat)[[1]][2]
      }
      get_unit_index_from_statename_v <- Vectorize(get_unit_index_from_statename)
      # convert to long format and output
      to_gather <- colnames(sims)[3:length(colnames(sims))][!c(colnames(sims)[3:length(colnames(sims))]%in%object@shared_covarnames)] # all columns except time and .id
      to_select <- c(colnames(sims)[1:2], "unit", "stateobs", "val")
      to_arrange <- c(colnames(sims)[1], "unit", "stateobs")
      gathered <- sims %>%
        tidyr::gather_(key="stateobs", val="val", to_gather) %>%
        dplyr::mutate(ui = get_unit_index_from_statename_v(stateobs))%>%
        dplyr::mutate(unit = unit_names(object)[as.integer(ui)]) %>%
        dplyr::select(to_select) %>%
        dplyr::arrange_(.dots = to_arrange)
      stateobstype <- sapply(gathered$stateobs,FUN=function(x) stringr::str_extract(x,unit_stateobs_pat))
      gathered$stateobstype <- stateobstype
      gathered <- gathered %>%
        dplyr::select(-stateobs) %>%
        tidyr::spread(key = stateobstype, value = val) %>%
        dplyr::rename(unitname = unit)
      return(gathered)
    }
    if(format=="spatPomps"){
      # add back spatPomp components into a list of spatPomps
      if(nsim > 1){
        sp.list <- vector(mode="list", length = nsim)
        for(i in 1:length(sims)){
          sp <- new("spatPomp",sims[[i]],
                    unit_covarnames = object@unit_covarnames,
                    shared_covarnames = object@shared_covarnames,
                    runit_measure = object@runit_measure,
                    dunit_measure = object@dunit_measure,
                    eunit_measure = object@eunit_measure,
                    munit_measure = object@munit_measure,
                    vunit_measure = object@vunit_measure,
                    unit_names=object@unit_names,
                    unitname=object@unitname,
                    unit_statenames=object@unit_statenames,
                    unit_obsnames = object@unit_obsnames)
          sp.list[[i]] <- sp
        }
        return(sp.list)
      } else{
        sp <- new("spatPomp",sims,
                  unit_covarnames = object@unit_covarnames,
                  shared_covarnames = object@shared_covarnames,
                  dunit_measure = object@dunit_measure,
                  runit_measure = object@runit_measure,
                  eunit_measure = object@eunit_measure,
                  munit_measure = object@munit_measure,
                  vunit_measure = object@vunit_measure,
                  unit_names=object@unit_names,
                  unitname=object@unitname,
                  unit_statenames=object@unit_statenames,
                  unit_obsnames = object@unit_obsnames)
        return(sp)
      }
    }
  }
)

##' @export
setMethod(
  "show",
  signature=signature(object="spatPomp"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)

##' @name logLik-girfd_spatPomp
##' @title logLik
##' @aliases logLik,girfd_spatPomp-method
##' @rdname logLik
##' @export
setMethod(
  "logLik",
  signature=signature(object="girfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-bpfilterd_spatPomp
##' @title logLik
##' @aliases logLik,bpfilterd_spatPomp-method
##' @rdname logLik
##' @export
setMethod(
  "logLik",
  signature=signature(object="bpfilterd_spatPomp"),
  definition=function(object)object@loglik
)


##' @name logLik-abfd.spatPomp
##' @title loglik
##' @aliases logLik,abfd.spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="abfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-abfird_spatPomp
##' @title loglik
##' @aliases logLik,abfird_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="abfird_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-igirfd_spatPomp
##' @title loglik
##' @aliases logLik,igirfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="igirfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name spatPomp_Csnippet
##' @title spatPomp_Csnippet
##' @rdname spatPomp_Csnippet
##' @export
setMethod(
  "spatPomp_Csnippet",
  signature=signature(object="character"),
  definition=function(object, unit_statenames,...){
    if(missing(unit_statenames))
      return(pomp::Csnippet(object))
    else{
      sn_inits_lhs <- paste("double *",unit_statenames, sep = "")
      sn_inits_rhs <- paste("&", unit_statenames,"1;",sep="")
      sn_inits_vec <- paste(sn_inits_lhs, sn_inits_rhs, sep = " = ")
      sn_inits <- paste0(sn_inits_vec, collapse = "\n")
    }
    all_inits <- paste(sn_inits, sep = "\n")
    full_csnippet <- paste(all_inits, object, sep = "\n")
    return(pomp::Csnippet(full_csnippet))
  }
)

