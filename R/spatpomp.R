spatpomp <- function (data, units, unit_index, times, covar, tcovar, ..., unit_dmeasure,
  unit_statenames, global_statenames, verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("spatpomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (missing(units))
    stop(ep,sQuote("units")," is a required argument",call.=FALSE)

  if (missing(times))
    stop(ep,sQuote("times")," is a required argument",call.=FALSE)

  if (missing(unit_dmeasure)) unit_dmeasure <- function(y,x,t,params,log=FALSE,d,...)
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
    upos_name <- names(data)[upos]


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
    tpos_name <- names(data)[tpos]

    # make data into unit by var by time
    temp_data <- abind::abind(split(data, data[tpos], drop = TRUE), along = 3)
    rownames(temp_data) <- temp_data[,upos,1]
    dat <- temp_data[,-c(upos,tpos),,drop = FALSE] # unit by var by time
    storage.mode(dat) <- "double" # TO DO: do i need to make the time axis 1:length(times) instead of actual times?

    # make covariates into unit by var by time (assuming user provides a df with columns time, unit, covariate1, covariate2,...)
    if(!missing(covar)){
      upos_cov <- match(upos_name, names(covar))
      tpos_cov <- match(tpos_name, names(covar))
      temp_cov <- abind::abind(split(covar, covar[tpos_cov], drop = TRUE), along = 3)
      rownames(temp_cov) <- temp_cov[,upos_cov,1]
      cov <- temp_cov[,-c(upos_cov,tpos_cov),,drop = FALSE] # unit by var by time
      storage.mode(cov) <- "double"
    } else{
      cov <- array(dim = c(0,0,0))
    }

    # get all combinations of unit statenames and units. Concatenate global statenames
    statenames <- c(do.call(paste0,expand.grid(unit_statenames, 1:length(units))), global_statenames)

    # units slot contains unique units. unit_index is an "ordering" of units
    units <- unique(data[[upos]])
    if(missing(unit_index)){
      unit_index <- units
      names(unit_index) <- 1:length(units)
    }

    po <- pomp(subset(data, data[[upos]]==data[[upos]][1], select = -c(upos)),
               times=times,
               dmeasure=unit_dmeasure,
               statenames=statenames,
               obsnames=colnames(dat),
               covar = subset(covar, covar[[upos_cov]] == covar[[upos_cov]][1], select = -c(upos_cov)),
               tcovar = tcovar,
               ...,
               verbose=verbose
               )

    spo <- new("spatpomp",po,unit_dmeasure=unit_dmeasure,units=units,unit_index=unit_index,
        unit_statenames=unit_statenames,global_statenames=global_statenames)
    spo@data <- dat
    spo@covar <- cov
    spo
  } else {
      stop(ep,sQuote("data"), " must be a data frame", call.=FALSE)
  }


}

