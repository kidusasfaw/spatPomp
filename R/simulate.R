##' Simulation of a spatiotemporal partially-observed Markov process
##'
##' \code{simulate} generates simulations of the latent and measurement
##' processes.
##'
##' @name simulate
##' @rdname simulate
##' @include spatPomp_class.R spatPomp.R
##'
##' @inheritParams pomp::simulate
##'
##' @param nsim number of simulations.
##' @param format the format of the simulated results. If the argument is
##' set to \code{'spatPomps'}, the default behavior, then the output is a
##' \code{list} of \code{spatPomp} objects. Options are \code{'spatPomps'}
##' and \code{'data.frame'}.
##' @importFrom utils data
##' @examples
##' # Complete examples are provided in the package tests
##' \dontrun{
##' # Get a spatPomp object
##' b <- bm(U=2, N=5)
##' # Get 2 simulations from same model as data.frame
##' sims <- simulate(b, nsim=2, format='data.frame')
##' }
NULL

setGeneric(
  "simulate",
  function (object, nsim=1, seed=NULL, ...)
    standardGeneric("simulate")
)

##' @name simulate-spatPomp
##' @aliases simulate,spatPomp-method
##' @rdname simulate
##' @return if \code{format='spatPomps'} and \code{nsim=1} an object of class \sQuote{spatPomp} representing a simulation from the model in \code{object} is returned.
##' If \code{format='spatPomps'} and \code{nsim>1} a list of class \sQuote{spatPomp} objects is returned.
##' If \code{format='data.frame'} then a class \sQuote{data.frame} object is returned.
##' @export
setMethod(
  "simulate",
  signature=signature(object="spatPomp"),
  definition=function(object, nsim = 1, seed = NULL,
    format = c("spatPomps", "data.frame"),
    include.data = FALSE,...) {
    format <- match.arg(format)
    if(format == 'spatPomps') sims <- pomp::simulate(pomp(object), format = 'pomps', nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format == 'data.frame') sims <- pomp::simulate(pomp(object), format = format, nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format=="data.frame"){
      get_unit_index_from_statename <- function(statename){
        stringr::str_extract(statename, "[[:digit:]]+$")
      }
      get_state_obs_type_from_statename <- function(statename) {
        gsub("[[:digit:]]+$", "", statename)
      }
      # convert to long format and output
      to_gather <- colnames(sims)[3:length(colnames(sims))][!c(colnames(sims)[3:length(colnames(sims))]%in%object@shared_covarnames)] # all columns except time and .id
      to_select <- c(colnames(sims)[1:2], "unitname", "stateobs", "val")
      to_arrange <- c(colnames(sims)[1], "unitname", "stateobs")
      sims %>%
        tidyr::pivot_longer(cols = to_gather, names_to = "stateobs",
	  values_to = "val") -> tmp
      tmp$ui <- get_unit_index_from_statename(tmp$stateobs)
      tmp$unitname <- unit_names(object)[as.integer(tmp$ui)]
      dplyr::select(tmp,dplyr::all_of(to_select)) %>%
        dplyr::arrange(dplyr::across(to_arrange)) -> tmp
      tmp$stateobstype <- get_state_obs_type_from_statename(tmp$stateobs)
      dplyr::select(tmp, -"stateobs") %>%
        tidyr::pivot_wider(names_from = "stateobstype",
	  values_from = 'val') -> gathered
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
