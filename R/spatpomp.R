spatpomp <- function (data, units, times, ..., unit_dmeasure,
  unit_statenames, global_statenames, verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("spatpomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (missing(units))
    stop(ep,sQuote("units")," is a required argument",call.=FALSE)

  if (missing(times))
    stop(ep,sQuote("times")," is a required argument",call.=FALSE)

  if (missing(unit_dmeasure)) unit_dmeasure <- function(y,x,times,params,d,log=FALSE,...)
    stop(sQuote("unit_dmeasure")," not specified")

  if (missing(unit_statenames)) unit_statenames <- character(0)

  if (missing(global_statenames)) global_statenames <- character(0)


  if (is.data.frame(data)) {
    ## 'data' is a data frame. find the units position
    if ((is.numeric(units) && (units<1 || units>ncol(data) ||
        units!=as.integer(units))) ||
        (is.character(units) && (!(units%in%names(data)))) ||
        (!is.numeric(units) && !is.character(units)) ||
        length(units)!=1) {
      stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("units"),
        " must identify a single column of ",sQuote("data"),
        " either by name or by index.",call.=FALSE)
    }
    if (is.numeric(units)) {
      upos <- as.integer(units)
    } else if (is.character(units)) {
      upos <- match(units,names(data))
    }

    if ((is.numeric(times) && (times<1 || times>ncol(data) ||times!=as.integer(times))) ||
        (is.character(times) && (!(times%in%names(data)))) ||
        (!is.numeric(times) && !is.character(times)) ||
        length(times)!=1) {
      stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("times"),
        " must identify a single column of ",sQuote("data"),
        " either by name or by index.",call.=FALSE)
    }
    if (is.numeric(times)) {
      tpos <- as.integer(times)
    } else if (is.character(times)) {
      tpos <- match(times,names(data))
    }
    units <- data[[upos]]
    temp_data <- abind::abind(split(data, data[tpos], drop = TRUE), along = 3)
    rownames(temp_data) <- temp_data[,upos,1]
    dat <- temp_data[,-c(upos,tpos),]
    storage.mode(dat) <- "numeric" # TO DO: do i need to make the time axis 1:length(times) instead of actual times?

    po <- pomp(subset(data, data[[upos]]==data[[upos]][1], select = -c(upos)),times=times,
      dmeasure=NULL,...,verbose=verbose)
    spo <- new("spatpomp",po,unit_dmeasure=unit_dmeasure,units=units,
        unit_statenames=unit_statenames,global_statenames=global_statenames)
    spo@data <- dat
    spo

  } else {
    stop(ep,sQuote("data"),
      " must be a data frame or an object of class ",sQuote("spatpomp"),
      call.=FALSE)
  }
}

