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
##' # Get a spatPomp object
##' b <- bm(U=5, N=10)
##' # Get 10 simulations from same model as data.frame
##' sims <- simulate(b, nsim=10, format='data.frame')
NULL

setGeneric(
  "simulate",
  function (object, nsim=1, seed=NULL, ...)
    standardGeneric("simulate")
)

##' @name simulate-spatPomp
##' @aliases simulate,spatPomp-method
##' @rdname simulate
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
        dplyr::mutate(ui = get_unit_index_from_statename_v(.data$stateobs))%>%
        dplyr::mutate(unit = unit_names(object)[as.integer(.data$ui)]) %>%
        dplyr::select(to_select) %>%
        dplyr::arrange_(.dots = to_arrange)
      stateobstype <- sapply(gathered$stateobs,FUN=function(x) stringr::str_extract(x,unit_stateobs_pat))
      gathered$stateobstype <- stateobstype
      gathered <- gathered %>%
        dplyr::select(-.data$stateobs) %>%
        tidyr::spread(key = stateobstype, value = .data$val) %>%
        dplyr::rename(unitname = .data$unit)
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
