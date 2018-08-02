#' @include spatpomp_class.R

spatpomp <- function (data, units, unit_index, times, covar, tcovar, t0, ...,
  unit_dmeasure, unit_statenames, global_statenames, rprocess, rmeasure,
  dprocess, dmeasure, initializer, cdir,cfile, shlib.args, userdata, PACKAGE,
  globals, statenames, paramnames, obstypes, zeronames, covarnames,
  verbose = getOption("verbose",FALSE)) {

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

    # get observation types
    obstypes <- names(data)[-c(upos,tpos)]

    # units slot contains unique units. unit_index is an "ordering" of units
    units <- unique(data[[upos]])
    if(missing(unit_index)){
      unit_index <- units
      names(unit_index) <- 1:length(units)
    }

    # make data into unit by var by time
    # temp_data <- abind::abind(split(data, data[tpos], drop = TRUE), along = 3)
    # rownames(temp_data) <- temp_data[,upos,1]
    # dat <- temp_data[,-c(upos,tpos),,drop = FALSE] # unit by var by time
    # storage.mode(dat) <- "double" # TO DO: do i need to make the time axis 1:length(times) instead of actual times?

    # make data into a dataframe that pomp would expect
    tmp <- names(unit_index)
    names(tmp) <- unit_index
    pomp_data <- data %>% mutate(ui = tmp[match(data[,upos_name], names(tmp))])
    pomp_data <- pomp_data %>% tidyr::gather(obstypes, key = 'obsname', value = 'val') %>% arrange(pomp_data[,tpos_name], obsname, ui)
    pomp_data <- pomp_data %>% dplyr::mutate(obsname = paste0(obsname,ui)) %>% dplyr::select(-upos) %>% dplyr::select(-ui)
    pomp_data <- pomp_data %>% tidyr::spread(key = obsname, value = val)


    # make covariates into unit by var by time (assuming user provides a df with columns time, unit, covariate1, covariate2,...)
    # if(!missing(covar)){
    #   upos_cov <- match(upos_name, names(covar))
    #   tpos_cov <- match(tpos_name, names(covar))
    #   temp_cov <- abind::abind(split(covar, covar[tpos_cov], drop = TRUE), along = 3)
    #   rownames(temp_cov) <- temp_cov[,upos_cov,1]
    #   cov <- temp_cov[,-c(upos_cov,tpos_cov),,drop = FALSE] # unit by var by time
    #   storage.mode(cov) <- "double"
    # } else{
    #   cov <- array(dim = c(0,0,0))
    # }

    # make covariates into a dataframe that pomp would expect
    if(!missing(covar)){
      upos_cov <- match(upos_name, names(covar))
      tpos_cov <- match(tpos_name, names(covar))
      covariate_names <- names(covar)[-c(upos_cov, tpos_cov)]
      tmp <- names(unit_index)
      names(tmp) <- unit_index
      pomp_covar <- covar %>% mutate(ui = match(covar[,upos_name], names(tmp)))
      pomp_covar <- pomp_covar %>% tidyr::gather(covariate_names, key = 'covname', value = 'val')
      pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(covname,ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-ui)
      pomp_covar <- pomp_covar %>% tidyr::spread(key = covname, value = val)
    }

    # get all combinations of unit statenames and units. Concatenate global statenames
    if(!missing(unit_statenames)){
      if(!missing(global_statenames)) pomp_statenames <- c(paste0(rep(unit_statenames,each=length(units)),1:length(units)),global_statenames)
      else pomp_statenames <- paste0(rep(unit_statenames,each=length(units)),1:length(units))
    } else pomp_statenames <- character(0)

    # get the observation names of the pomp dataframe.
    obsnames <- names(pomp_data)[-c(tpos)]
    if(missing(globals)) globals <- Csnippet(paste0("int nunits = ",length(units),";\n"))
    else globals <- Csnippet(paste0(paste0("\nint nunits = ",length(units),";\n"),globals@text))

    # create the pomp object
    po <- pomp(data = pomp_data,
      rprocess = rprocess,
      rmeasure = rmeasure,
      dprocess = dprocess,
      dmeasure = dmeasure,
      initializer = initializer,
      times=times,
      statenames=pomp_statenames,
      zeronames=zeronames,
      covar = pomp_covar,
      tcovar = tcovar,
      paramnames = paramnames,
      globals = globals,
      t0 = t0,
      ...,
      verbose=verbose
    )
    ## preliminary error checking
    if (missing(cdir)) cdir <- NULL
    if (missing(cfile)) cfile <- NULL
    if (missing(shlib.args)) shlib.args <- NULL
    if (missing(userdata)) userdata <- list()
    added.userdata <- list(...)
    if (length(added.userdata)>0) {
      message("In ",sQuote("spatpomp"),
        ": the following unrecognized argument(s) ",
        "will be stored for use by user-defined functions: ",
        paste(sQuote(names(added.userdata)),collapse=","))
      userdata[names(added.userdata)] <- added.userdata
    }

    ## name of shared object library
    if (missing(PACKAGE)) PACKAGE <- NULL
    PACKAGE <- as.character(PACKAGE)

    if (missing(globals)) globals <- NULL
    if (!is(globals,"Csnippet"))
      globals <- as.character(globals)

    ## defaults for names of states, parameters, and accumulator variables
    statenames <- pomp_statenames
    if (missing(paramnames)) paramnames <- character(0)
    if (missing(zeronames)) zeronames <- character(0)
    statenames <- as.character(statenames)
    paramnames <- as.character(paramnames)
    zeronames <- as.character(zeronames)

    ## check for duplicate names
    if (anyDuplicated(statenames)) {
      stop("all ",sQuote("statenames")," must be unique", call.=FALSE)
    }
    if (anyDuplicated(paramnames)) {
      stop("all ",sQuote("paramnames")," must be unique", call.=FALSE)
    }
    if (anyDuplicated(zeronames)) {
      stop("all ",sQuote("zeronames")," must be unique", call.=FALSE)
    }

    # arrange covariates
    if (missing(covarnames) || length(covarnames)==0)
      covarnames <- as.character(colnames(pomp_covar))
    if (!all(covarnames%in%colnames(pomp_covar))) {
      missing <- covarnames[!(covarnames%in%colnames(covar))]
      stop("covariate(s) ",
        paste(sapply(missing,sQuote),collapse=","),
        " are not among the columns of ",sQuote("covar"),call.=FALSE)
    }
    ## handle unit_dmeasure C Snippet
    ud_template <- list(
      unit_dmeasure=list(
        slotname="unit_dmeasure",
        Cname="__spatpomp_unit_dmeasure",
        proto=quote(unit_dmeasure(y,x,t,d,params,log,...)),
        header="\nvoid __spatpomp_unit_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
        footer="\n}\n\n",
        vars=list(
          params=list(
            names=quote(paramnames),
            cref="__p[__parindex[{%v%}]]"
          ),
          covars=list(
            names=quote(covarnames),
            cref="__covars[__covindex[{%v%}]]"
          ),
          states=list(
            names=quote(statenames),
            cref="__x[__stateindex[{%v%}*nunits+unit-1]]"
          ),
          obs=list(
            names=quote(obsnames),
            cref="__y[__obsindex[{%v%}*nunits+unit-1]]"
          ),
          lik=list(
            names="lik",
            cref="__lik[0]"
          )
        )
      )
    )

    hitches <- pomp::hitch(
      unit_dmeasure=unit_dmeasure,
      templates=ud_template,
      obsnames=obstypes,
      statenames=unit_statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      PACKAGE=PACKAGE,
      globals=globals,
      cfile=cfile,
      cdir=cdir,
      shlib.args=shlib.args,
      verbose=verbose
    )

    pomp:::solibs(po) <- hitches$lib

    new("spatpomp",po,
      unit_dmeasure=hitches$funs$unit_dmeasure,
      units=units,
      unit_index=unit_index,
      unit_statenames=unit_statenames,
      global_statenames=global_statenames,
      obstypes = obstypes)

  } else {
    stop(ep,sQuote("data"), " must be a data frame", call.=FALSE)
  }
}
