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
        unit_dmeasure=from@dmeasure,
        units="unit",
        unit_statenames=character(0),
        obstypes = rownames(from@data))
  }
)

##' @export
setMethod(
  "spat_units",
  signature=signature(x="spatPomp"),
  definition=function(x,...) x@units
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
  "plot",
  signature=signature(x="spatPomp"),
  definition=function (x, ...) {
    plot(as.data.frame(x))
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
      # convert to long format and output
      to.gather <- colnames(sims)[3:length(colnames(sims))] # all columns except time and .id
      to.select <- c(colnames(sims)[1:2], "unit", "stateobs", "val")
      to.arrange <- c(colnames(sims)[1], "unit", "stateobs")
      gathered <- sims %>% tidyr::gather_(key="stateobs", val="val", to.gather) %>%
        dplyr::mutate(ui = stringr::str_extract(stateobs,"[0-9]+"))%>%
        dplyr::mutate(unit = spat_units(object)[ui])%>%
        dplyr::select_(.dots = to.select) %>%
        dplyr::arrange_(.dots = to.arrange)
      stateobstype <- sapply(gathered$stateobs,FUN=function(x) stringr::str_split(x,"[0-9]+")[[1]][1])
      gathered$stateobstype <- stateobstype
      gathered <- gathered %>%
        dplyr::select(-stateobs) %>%
        tidyr::spread(key = stateobstype, value = val)
      return(gathered)
    }
    if(format=="spatPomps"){
      # add back spatPomp components into a list of spatPomps
      if(nsim > 1){
        sp.list <- vector(mode="list", length = nsim)
        for(i in 1:length(sims)){
          sp <- new("spatPomp",sims[[i]],
                    unit_rmeasure = object@unit_rmeasure,
                    unit_dmeasure = object@unit_dmeasure,
                    unit_emeasure = object@unit_emeasure,
                    unit_mmeasure = object@unit_mmeasure,
                    unit_vmeasure = object@unit_vmeasure,
                    units=object@units,
                    unit_statenames=object@unit_statenames,
                    obstypes = object@obstypes)
          sp.list[[i]] <- sp
        }
        return(sp.list)
      } else{
        sp <- new("spatPomp",sims,
                  unit_dmeasure = object@unit_dmeasure,
                  unit_rmeasure = object@unit_rmeasure,
                  unit_emeasure = object@unit_emeasure,
                  unit_mmeasure = object@unit_mmeasure,
                  unit_vmeasure = object@unit_vmeasure,
                  units=object@units,
                  unit_statenames=object@unit_statenames,
                  obstypes = object@obstypes)
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
##' @title loglik
##' @aliases logLik,girfd_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="girfd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-asifd.spatPomp
##' @title loglik
##' @aliases logLik,asifd.spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="asifd_spatPomp"),
  definition=function(object)object@loglik
)

##' @name logLik-asifird_spatPomp
##' @title loglik
##' @aliases logLik,asifird_spatPomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="asifird_spatPomp"),
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
  definition=function(object, unit_statenames, unit_covarnames){
    if(missing(unit_statenames) && missing(unit_covarnames))
      return(pomp::Csnippet(object))
    else{
      if(missing(unit_statenames)) sn.inits <- character()
      else{
        sn.inits.lhs <- paste("double *",unit_statenames, sep = "")
        sn.inits.rhs <- paste("&", unit_statenames,"1;",sep="")
        sn.inits.vec <- paste(sn.inits.lhs, sn.inits.rhs, sep = " = ")
        sn.inits <- paste0(sn.inits.vec, collapse = "\n")
      }
      if(missing(unit_covarnames)) cn.inits <- character()
      else{
        cn.inits.lhs <- paste("const double *",unit_covarnames, sep = "")
        cn.inits.rhs <- paste("&", unit_covarnames,"1;",sep="")
        cn.inits.vec <- paste(cn.inits.lhs, cn.inits.rhs, sep = " = ")
        cn.inits <- paste0(cn.inits.vec, collapse = "\n")
      }
      all.inits <- paste(sn.inits, cn.inits, sep = "\n")
      full.csnippet <- paste(all.inits, object, sep = "\n")
      return(pomp::Csnippet(full.csnippet))
    }
  }
)

