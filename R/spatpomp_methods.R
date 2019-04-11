#' @include spatpomp_class.R
#'
## this file contains some basic methods definitions
setAs(
  from="pomp",
  to="spatpomp",
  def = function (from) {
    new("spatpomp",from,
        unit_dmeasure=from@dmeasure,
        units="unit",
        unit_index="unit",
        unit_statenames=character(0),
        global_statenames=character(0),
        obstypes = rownames(from@data))
  }
)
## extract the vector of units
setMethod(
  "unit",
  signature=signature(x="spatpomp"),
  definition=function(x,...) x@units
)

## extract the unit index
setMethod(
  "unit_ix",
  signature=signature(x="spatpomp"),
  definition=function(x,...) x@unit_index
)


setMethod(
  "print",
  signature=signature(x="spatpomp"),
  definition=function (x, ...) {
    cat("<object of class ",sQuote("spatpomp"),">\n",sep="")
    invisible(x)
  }
)

setMethod(
  "simulate",
  signature=signature(object="spatpomp"),
  definition=function(object, nsim = 1, seed = NULL,
                       format = c("spatpomps", "arrays", "data.frame"),
                       include.data = FALSE,...) {
    format <- match.arg(format)
    if(format == 'spatpomps') sims <- simulate(pomp(object), format = 'pomps', nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format == 'data.frame') sims <- simulate(pomp(object), format = format, nsim = nsim, include.data = include.data, seed = seed, ...)
    if(format=="data.frame"){
      # convert to long format and output
      to.gather <- colnames(sims)[3:length(colnames(sims))] # all columns except time and .id
      to.select <- c(colnames(sims)[1:2], "unit", "stateobs", "val")
      to.arrange <- c(colnames(sims)[1], "unit", "stateobs")
      gathered <- sims %>% tidyr::gather_(key="stateobs", val="val", to.gather) %>%
        dplyr::mutate(ui = stringr::str_extract(stateobs,"[0-9]+"))%>%
        dplyr::mutate(unit = object@unit_index[ui])%>%
        dplyr::select_(.dots = to.select) %>%
        dplyr::arrange_(.dots = to.arrange)
      stateobstype <- sapply(gathered$stateobs,FUN=function(x) stringr::str_split(x,"[0-9]+")[[1]][1])
      gathered$stateobstype <- stateobstype
      gathered <- gathered %>%
        dplyr::select(-stateobs) %>%
        tidyr::spread(key = stateobstype, value = val)
      return(gathered)
    }
    if(format=="spatpomps"){
      # add back spatpomp components into a list of spatpomps
      if(nsim > 1){
        sp.list <- vector(mode="list", length = nsim)
        for(i in 1:length(sims)){
          sp <- new("spatpomp",sims[[i]],
                    unit_rmeasure = object@unit_rmeasure,
                    unit_dmeasure = object@unit_dmeasure,
                    units=object@units,
                    unit_index=object@unit_index,
                    unit_statenames=object@unit_statenames,
                    global_statenames=object@global_statenames,
                    obstypes = object@obstypes)
          sp.list[[i]] <- sp
        }
        return(sp.list)
      } else{
        sp <- new("spatpomp",sims,
                  unit_dmeasure = object@unit_dmeasure,
                  units=object@units,
                  unit_index=object@unit_index,
                  unit_statenames=object@unit_statenames,
                  global_statenames=object@global_statenames,
                  obstypes = object@obstypes)
        return(sp)
      }
    }
  }
)

setMethod(
  "show",
  signature=signature(object="spatpomp"),
  definition=function (object) {
    cat("<object of class ",sQuote(as.character(class(object))),">\n",sep="")
    invisible(NULL)
  }
)

