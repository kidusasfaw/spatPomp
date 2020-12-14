##' Basic methods for the spatPomp class
##'
##' @name spatpomp-methods
##' @title spatpomp-methods
##' @rdname spatPomp-methods
##'
##' @include spatPomp_class.R
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

setGeneric("unit_names", function(x,...)standardGeneric("unit_names"))

##' Unit names of a spatiotemporal model
##'
##' \code{unit_names} outputs the contents of the \code{unit_names} slot
##' of a \code{spatPomp} object. The order in which the spatial units
##' appear in the output vector determines the order in which latent
##' states and observations for the spatial units are stored.
##'
##' @name unit_names
##' @rdname unit_names
##' @include spatPomp_class.R spatPomp.R
##'
##' @inheritParams abf
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
##' Simulation of a spatiotemporal partially-observed Markov process
##'
##' \code{simulate} generates simulations of the latent and measurement
##' processes.
##'
##' @name simulate
##' @rdname simulate
##' @include spatPomp_class.R spatPomp.R
##'
##' @inheritParams abf
##'
##' @param nsim number of simulations.
##' @param format the format of the simulated results. If the argument is
##' set to \code{'spatPomps'}, the default behavior, then the output is a
##' \code{list} of \code{spatPomp} objects. Options are \code{'spatPomps'}
##' and \code{'data.frame'}.
##' @examples
##' # Get a spatPomp object
##' b <- bm(U=5, N=10)
##' # Get 10 simulations from same model as data.frame
##' sims <- simulate(b, nsim=10, format='data.frame')
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
                    unit_accumvars = object@unit_accumvars,
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
                  unit_accumvars = object@unit_accumvars,
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


##' @name logLik-abfd_spatPomp
##' @title loglik
##' @aliases logLik,abfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="abfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-iabfd_spatPomp
##' @title loglik
##' @aliases logLik,iabfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="iabfd_spatPomp"),
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

setGeneric("spatPomp_Csnippet", function(code,...)standardGeneric("spatPomp_Csnippet"))

##' C snippets
##'
##' \code{spatPomp_Csnippet()} is used to provide snippets of \proglang{C}
##' code that specify model components. It functions similarly to \code{Csnippet()} from
##' the \pkg{pomp} package; in fact, the output of \code{spatPomp_Csnippet} is an object
##' of class \code{Csnippet}.  It additionally provides some arguments that allow the user
##' to stay focused on model development in the spatiotemporal context  where
##' model size grows.
##'
##' @name spatPomp_Csnippet
##' @rdname spatPomp_Csnippet
##' @include spatPomp_class.R spatPomp.R
##'
##' @param code encodes a component of a spatiotemporal POMP model using \proglang{C} code
##' @param unit_statenames a subset of the \code{unit_statenames} slot of
##' the \code{spatPomp} object for which we are writing a model. This argument
##' allows the user to get variables that can be indexed conveniently to update
##' states and measurements in a loop. See examples for more details.
##' @param unit_obsnames a subset of the \code{unit_obsnames} slot of
##' the \code{spatPomp} object for which we are writing a model. This argument
##' allows the user to get variables that can be indexed conveniently to update
##' states and measurements in a loop. See examples for more details.
##' @param unit_covarnames if the model has covariate information for each unit,
##' the names of the covariates for each unit can be supplied to this argument.
##' This allows the user to get variables that can be indexed conveniently to
##' use incorporate the covariate information in a loop. See examples for more
##' details.
##' @param unit_ivpnames This argument is particularly useful when specifying the
##' \code{rinit} model component. The \code{paramnames} argument to the
##' \code{spatPomp()} constructor often has names for initial value
##' parameters for the latent states (e.g. \code{S1_0}, \code{S2_0} for the
##' the quantity of susceptibles at unit 1 and unit 2 at the initial time in an
##' SIR model). By supplying \code{unit_ivpnames}, we can get variables
##' that can be easily indexed to reference the initial value parameters (in
##' the previous example, \code{unit_ivpnames=c('S')} we can get a variable
##' named \code{S_0} that we can index as \code{S_0[0]} and \code{S_0[1]} to
##' refer to \code{S1_0} and \code{S2_0}). See examples for more details.
##'
##' @param unit_vfnames This argument is particularly useful when specifying the
##' \code{skeleton} model component. For all components of the latent state,
##' the user can assume a variable defining the time-derivative is pre-defined (e.g.
##' \code{DS1} and \code{DS2} for the time-derivative of the quantity of the
##' susceptibles at unit 1 and unit 2 in an SIR model). By supplying
##' \code{unit_vfnames}, we can get variables that can be easily indexed to
##' reference these variables (in the previous example,
##' setting \code{unit_vfnames=c('S')} gets us a variable
##' named \code{DS} that we can index as \code{DS[0]} and \code{DS[1]} to
##' refer to \code{DS1} and \code{DS2}). See examples for more details.
##'
##' @examples
##' # Set initial states for Brownian motion
##' bm_rinit <- spatPomp_Csnippet(
##'   unit_statenames = c("X"),
##'   unit_ivpnames = c("X"),
##'   code = "
##'     for (int u = 0; u < U; u++) {
##'       X[u]=X_0[u];
##'     }
##'   "
##' )
##' # Skeleton for Brownian motion
##' bm_skel <- spatPomp_Csnippet(
##'   unit_statenames = c("X"),
##'   unit_vfnames = c("X"),
##'   code = "
##'       for (int u = 0 ; u < U ; u++) {
##'         DX[u] = 0;
##'       }
##'   "
##')
##' @export
setMethod(
  "spatPomp_Csnippet",
  signature=signature(code="character"),
  definition=function(code, unit_statenames, unit_obsnames, unit_covarnames, unit_ivpnames, unit_paramnames, unit_vfnames){
    if(missing(unit_statenames) &&
       missing(unit_obsnames) &&
       missing(unit_covarnames) &&
       missing(unit_ivpnames) &&
       missing(unit_paramnames) &&
       missing(unit_vfnames))
      return(pomp::Csnippet(code))
    sn_inits <- on_inits <- cn_inits <- in_inits <- pn_inits <- vn_inits <- character()
    if(!missing(unit_statenames)){
      sn_inits_lhs <- paste("double *",unit_statenames, sep = "")
      sn_inits_rhs <- paste("&", unit_statenames,"1;",sep="")
      sn_inits_vec <- paste(sn_inits_lhs, sn_inits_rhs, sep = " = ")
      sn_inits <- paste0(sn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_obsnames)){
      on_inits_lhs <- paste("const double *",unit_obsnames, sep = "")
      on_inits_rhs <- paste("&", unit_obsnames,"1;",sep="")
      on_inits_vec <- paste(on_inits_lhs, on_inits_rhs, sep = " = ")
      on_inits <- paste0(on_inits_vec, collapse = "\n")
    }
    if(!missing(unit_covarnames)){
      cn_inits_lhs <- paste("const double *",unit_covarnames, sep = "")
      cn_inits_rhs <- paste("&", unit_covarnames,"1;",sep="")
      cn_inits_vec <- paste(cn_inits_lhs, cn_inits_rhs, sep = " = ")
      cn_inits <- paste0(cn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_ivpnames)){
      left_ivpnames = paste(unit_ivpnames, "_0", sep = "")
      in_inits_lhs <- paste("const double *",left_ivpnames, sep = "")
      in_inits_rhs <- paste("&", unit_ivpnames,"1_0;",sep="")
      in_inits_vec <- paste(in_inits_lhs, in_inits_rhs, sep = " = ")
      in_inits <- paste0(in_inits_vec, collapse = "\n")
    }
    if(!missing(unit_paramnames)){
      pn_inits_lhs <- paste("const double *",unit_paramnames, sep = "")
      pn_inits_rhs <- paste("&", unit_paramnames,"1;",sep="")
      pn_inits_vec <- paste(pn_inits_lhs, pn_inits_rhs, sep = " = ")
      pn_inits <- paste0(pn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_vfnames)){
      vn_inits_lhs <- paste("double *D",unit_vfnames, sep = "")
      vn_inits_rhs <- paste("&", unit_vfnames,"1;",sep="")
      vn_inits_vec <- paste(vn_inits_lhs, vn_inits_rhs, sep = " = ")
      vn_inits <- paste0(vn_inits_vec, collapse = "\n")
    }
    all_inits <- paste(sn_inits, on_inits, cn_inits, in_inits, pn_inits, vn_inits, sep = "\n")
    full_csnippet <- paste(all_inits, code, sep = "\n")
    return(pomp::Csnippet(full_csnippet))
  }
)

