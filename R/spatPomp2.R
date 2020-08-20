##' @export
spatPomp2 <- function (data, units, times, covar, tcovar, t0, ...,
                      unit_emeasure, unit_mmeasure, unit_vmeasure, unit_dmeasure, unit_rmeasure,
                      unit_statenames, rprocess, rmeasure,
                      dprocess, dmeasure, skeleton, rinit, cdir,cfile, shlib.args, PACKAGE,
                      globals, statenames, paramnames, params, unit_obsnames, accumvars, covarnames,
                      partrans, verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("spatPomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (!inherits(data,what=c("data.frame","spatPomp")))
    pStop("pomp",sQuote("data")," must be a data frame or an object of ",
          "class ",sQuote("spatPomp"),".")

  ## return as quickly as possible if no work is to be done
  if (is(data,"spatPomp") && missing(times) && missing(t0) &&
      missing(unit_dmeasure) && missing(unit_emeasure) &&
      missing(unit_vmeasure) && missing(unit_mmeasure) &&
      missing(unit_rmeasure) &&
      missing(rinit) && missing(rprocess) && missing(dprocess) &&
      missing(rmeasure) && missing(dmeasure) && missing(skeleton) &&
      missing(rprior) && missing(dprior) && missing(partrans) &&
      missing(covar) && missing(tcovar) && missing(params) && missing(accumvars) &&
      length(list(...)) == 0)
    return(as(data,"spatPomp"))

  if (missing(times)) times <- NULL
  if (missing(units)) units <- NULL
  tryCatch(
    construct_spatPomp(
      data=data,times=times,units=units,t0=t0,...,
      rinit=rinit,rprocess=rprocess,dprocess=dprocess,
      rmeasure=rmeasure,dmeasure=dmeasure,
      skeleton=skeleton,rprior=rprior,dprior=dprior,partrans=partrans,
      params=params,covar=covar,tcovar=tcovar,accumvars=accumvars,
      unit_dmeasure=unit_dmeasure,unit_emeasure=unit_emeasure,
      unit_vmeasure=unit_vmeasure,unit_mmeasure=unit_mmeasure,
      unit_rmeasure=unit_rmeasure,unit_statenames=unit_statenames,
      obsnames=obsnames,statenames=statenames,paramnames=paramnames,
      covarnames=covarnames,PACKAGE=PACKAGE,
      globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args,
      compile=compile, verbose=verbose
    ),
    error = function (e) pomp:::pStop_(conditionMessage(e))
  )
}

setMethod(
  "construct_spatPomp",
  signature=signature(data="data.frame", times="character", units="character"),
  definition = function (data, times, units, t0, ...,
                         rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
                         partrans, params, covar, tcovar, accumvars, unit_dmeasure, unit_emeasure,
                         unit_vmeasure, unit_mmeasure, unit_rmeasure, unit_statenames, obsnames,
                         statenames, paramnames, covarnames, PACKAGE, globals,
                         cdir, cfile, shlib.args, compile, verbose) {

    if (anyDuplicated(names(data)))
      pStop_("names of data variables must be unique.")

    if (missing(t0)) reqd_arg(NULL,"t0")

    tpos <- match(times,names(data),nomatch=0L)
    upos <- match(units,names(data),nomatch=0L)

    if (length(times) != 1 || tpos == 0L)
      pomp:::pStop_(sQuote("times")," does not identify a single column of ",
             sQuote("data")," by name.")
    if (length(units) != 1 || upos == 0L)
      pomp:::pStop_(sQuote("units")," does not identify a single column of ",
             sQuote("data")," by name.")

    timename <- times
    unitname <- units

    # units slot contains unique units. unit_index is an "ordering" of units
    units <- unique(data[[upos]])
    unit_index <- units
    names(unit_index) <- 1:length(units)

    # get observation types
    unit_obsnames <- names(data)[-c(upos,tpos)]

    # if missing workhorses, set to default
    if (missing(rinit)) rinit <- NULL

    if (missing(rprocess) || is.null(rprocess)) {
      rprocess <- pomp:::rproc_plugin()
    }

    if (missing(dprocess)) dprocess <- NULL
    if (missing(rmeasure)) rmeasure <- NULL
    if (missing(dmeasure)) dmeasure <- NULL
    if (missing(unit_dmeasure)) unit_dmeasure <- NULL
    if (missing(unit_emeasure)) unit_emeasure <- NULL
    if (missing(unit_vmeasure)) unit_vmeasure <- NULL
    if (missing(unit_mmeasure)) unit_mmeasure <- NULL
    if (missing(unit_rmeasure)) unit_rmeasure <- NULL

    if (missing(skeleton) || is.null(skeleton)) {
      skeleton <- pomp:::skel_plugin()
    }

    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL

    if (missing(partrans) || is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(params)) params <- numeric(0)
    if (is.list(params)) params <- unlist(params)

    if (missing(covar)) covar <- pomp::covariate_table()

    # make data into a dataframe that pomp would expect
    tmp <- names(unit_index)
    names(tmp) <- unit_index
    pomp_data <- data %>% dplyr::mutate(ui = tmp[match(data[,unitname], names(tmp))])
    pomp_data <- pomp_data %>% tidyr::gather(unit_obsnames, key = 'obsname', value = 'val') %>% dplyr::arrange(pomp_data[,timename], obsname, ui)
    pomp_data <- pomp_data %>% dplyr::mutate(obsname = paste0(obsname,ui)) %>% dplyr::select(-upos) %>% dplyr::select(-ui)
    pomp_data <- pomp_data %>% tidyr::spread(key = obsname, value = val)
    dat_col_order <- vector(length = length(units)*length(unit_obsnames))
    for(ot in unit_obsnames){
      for(i in 1:length(units)){
        dat_col_order[i] = paste0(ot, i)
      }
    }
    pomp_data <- pomp_data[, c(timename, dat_col_order)]

    # make covariates into a dataframe that pomp would expect
    if(!missing(covar) && !missing(tcovar)){
      upos_cov <- match(unitname, names(covar))
      tpos_cov <- match(timename, names(covar))
      covariate_names <- names(covar)[-c(upos_cov, tpos_cov)]
      tmp <- names(unit_index)
      names(tmp) <- unit_index
      pomp_covar <- covar %>% dplyr::mutate(ui = match(covar[,unitname], names(tmp)))
      pomp_covar <- pomp_covar %>% tidyr::gather(covariate_names, key = 'covname', value = 'val')
      pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(covname,ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-ui)
      pomp_covar <- pomp_covar %>% tidyr::spread(key = covname, value = val)
      cov_col_order <- c()
      for(cn in covariate_names){
        for(i in 1:length(units)){
          cov_col_order = c(cov_col_order, paste0(cn, i))
        }
      }
      pomp_covar <- pomp_covar[, c(timename, cov_col_order)]
      # construct call of covariate_table function
      call_to_covar = list()
      call_to_covar[[1]] <- as.symbol("covariate_table")
      for(col in names(pomp_covar)){
        call_to_covar$col <- pomp_covar[,col]
      }
      call_to_covar$'times'=timename
    } else {
      pomp_covar <- NULL
      tcovar <- NULL
      call_to_covar <- NULL
    }
    if(is.null(tcovar)){
      cov.t <- NULL
    } else {
      cov.t <- covariate_table(pomp_covar, times = tcovar)
    }

    # get all combinations of unit statenames and units.
    if(!missing(unit_statenames)){
      pomp_statenames <- paste0(rep(unit_statenames,each=length(units)),1:length(units))
    # } else pomp_statenames <- NULL
    }

    # We will always have a global giving us the number of spatial units
    if(missing(globals)) globals <- Csnippet(paste0("const int nunits = ",length(units),";\n"))
    else globals <- Csnippet(paste0(paste0("\nconst int nunits = ",length(units),";\n"),globals@text))

    # create the pomp object
    po <- pomp(data = pomp_data,
               times=times,
               t0 = t0,
               rprocess = rprocess,
               rmeasure = rmeasure,
               dprocess = dprocess,
               dmeasure = dmeasure,
               skeleton = skeleton,
               rinit = rinit,
               statenames=pomp_statenames,
               accumvars=accumvars,
               covar = cov.t,
               paramnames = paramnames,
               globals = globals,
               cdir = cdir,
               cfile = cfile,
               shlib.args = shlib.args,
               partrans = partrans,
               ...,
               verbose=verbose
    )

    ## defaults for names of states, parameters, and accumulator variables
    if (missing(statenames)) statenames <- NULL
    if (missing(paramnames)) paramnames <- NULL
    if (missing(accumvars)) accumvars <- NULL
    accumvars <- unique(as.character(accumvars))
    if (!missing(paramnames)) mparamnames <- paste("M_", paramnames, sep = "")

    hitches <- pomp::hitch(
      unit_emeasure=unit_emeasure,
      unit_mmeasure=unit_mmeasure,
      unit_vmeasure=unit_vmeasure,
      unit_dmeasure=unit_dmeasure,
      unit_rmeasure=unit_rmeasure,
      templates=eval(spatPomp_workhorse_templates),
      obsnames = paste0(unit_obsnames,"1"),
      statenames = paste0(unit_statenames,"1"),
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
  new("spatPomp",po,
      unit_emeasure=hitches$funs$unit_emeasure,
      unit_mmeasure=hitches$funs$unit_mmeasure,
      unit_vmeasure=hitches$funs$unit_vmeasure,
      unit_dmeasure=hitches$funs$unit_dmeasure,
      unit_rmeasure=hitches$funs$unit_rmeasure,
      units=units,
      unit_statenames=unit_statenames,
      unit_obsnames = unit_obsnames)

  }
)

setMethod(
  "construct_spatPomp",
  signature=signature(data="spatPomp", times="NULL", units="NULL"),
  definition = function (data, times, units, t0, timename, ...,
                         rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
                         partrans, params, covar, tcovar, accumvars, unit_dmeasure, unit_emeasure,
                         unit_vmeasure, unit_mmeasure, unit_rmeasure, unit_statenames, obsnames,
                         statenames, paramnames, covarnames, PACKAGE, globals,
                         cdir, cfile, shlib.args, compile,.userdata, .solibs, verbose) {
    times <- data@times
    units <- data@unit_names
    unit_statenames <- data@unit_statenames
    unit_obsnames <- data@unit_obsnames
    if (missing(timename) || is.null(timename))
      timename <- "time"
    else
      timename <- as.character(timename)
    ## defaults for names of parameters, and accumulator variables
    if (!missing(unit_mmeasure) && missing(paramnames)){
      stop("Replacing  ",
           sQuote("unit_mmeasure"),
           " requires supplying ",sQuote("paramnames"),call.=FALSE)
    }
    if (missing(paramnames)) paramnames <- NULL
    if (missing(accumvars)) accumvars <- NULL
    accumvars <- unique(as.character(accumvars))
    if (!missing(paramnames)) mparamnames <- paste("M_", paramnames, sep = "")
    if (missing(covarnames)) covarnames <- NULL
    if (missing(.userdata)) .userdata <- list()
    added.userdata <- list(...)
    if (length(added.userdata)>0) {
      message("in ",sQuote("spatPomp"),": the unrecognized ",
              ngettext(length(added.userdata),"argument","arguments")," ",
              paste(sQuote(names(added.userdata)),collapse=","),
              ngettext(length(added.userdata)," is"," are"),
              " available for use by the POMP basic components."
      )
      .userdata[names(added.userdata)] <- added.userdata
    }
    # We will always have a global giving us the number of spatial units
    if(missing(globals)) globals <- Csnippet(paste0("const int nunits = ",length(units),";\n"))
    else globals <- Csnippet(paste0(globals@text,paste0("\nconst int nunits = ",length(units),";\n")))


    if (missing(t0)) t0 <- data@t0
    if (missing(rinit)) rinit <- data@rinit
    if (missing(rprocess)) rprocess <- data@rprocess
    else if (is.null(rprocess)) rprocess <- pomp:::rproc_plugin()
    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure
    if (missing(unit_dmeasure)) unit_dmeasure <- data@unit_dmeasure
    if (missing(unit_mmeasure)) unit_mmeasure <- data@unit_mmeasure
    if (missing(unit_vmeasure)) unit_vmeasure <- data@unit_vmeasure
    if (missing(unit_emeasure)) unit_emeasure <- data@unit_emeasure
    if (missing(unit_rmeasure)) unit_rmeasure <- data@unit_rmeasure
    if (missing(skeleton)) skeleton <- data@skeleton
    else if (is.null(skeleton)) skeleton <- skel_plugin()

    if (missing(rprior)) rprior <- data@rprior
    if (missing(dprior)) dprior <- data@dprior
    if (missing(partrans)) partrans <- data@partrans
    else if (is.null(partrans)) partrans <- parameter_trans()

    if (missing(params)) params <- data@params
    if (missing(covar)) covar <- data@covar
    if (missing(accumvars)) accumvars <- data@accumvars
    if (missing(.solibs)) .solibs <- data@solibs

    hitches <- pomp::hitch(
      rinit=rinit,
      step.fn=rprocess@step.fn,
      rate.fn=rprocess@rate.fn,
      dprocess=dprocess,
      rmeasure=rmeasure,
      dmeasure=dmeasure,
      rprior=rprior,
      dprior=dprior,
      toEst=partrans@to,
      fromEst=partrans@from,
      skeleton=skeleton@skel.fn,
      unit_emeasure=unit_emeasure,
      unit_mmeasure=unit_mmeasure,
      unit_vmeasure=unit_vmeasure,
      unit_dmeasure=unit_dmeasure,
      unit_rmeasure=unit_rmeasure,
      templates=eval(spatPomp_workhorse_templates),
      obsnames = paste0(unit_obsnames,"1"),
      statenames = paste0(unit_statenames,"1"),
      paramnames=paramnames,
      covarnames=covarnames,
      PACKAGE=PACKAGE,
      globals=globals,
      cfile=cfile,
      cdir=cdir,
      shlib.args=shlib.args,
      verbose=verbose
    )

    new(
      "spatPomp",
      data = data@data,
      times = times,
      units = units,
      unit_statenames = unit_statenames,
      unit_obsnames = unit_obsnames,
      t0 = t0,
      timename = timename,
      rinit = hitches$funs$rinit,
      rprocess = pomp:::rproc_plugin(
        rprocess,
        step.fn=hitches$funs$step.fn,
        rate.fn=hitches$funs$rate.fn
      ),
      dprocess = hitches$funs$dprocess,
      dmeasure = hitches$funs$dmeasure,
      rmeasure = hitches$funs$rmeasure,
      skeleton = pomp:::skel_plugin(
        skeleton,
        skel.fn=hitches$funs$skeleton
      ),
      dprior = hitches$funs$dprior,
      rprior = hitches$funs$rprior,
      partrans = pomp:::parameter_trans(
        toEst=hitches$funs$toEst,
        fromEst=hitches$funs$fromEst
      ),
      params = params,
      covar = covar,
      accumvars = accumvars,
      solibs = if (is.null(hitches$lib)) {
        .solibs
      } else {
        c(list(hitches$lib),.solibs)
      },
      userdata = .userdata,
      unit_emeasure=hitches$funs$unit_emeasure,
      unit_mmeasure=hitches$funs$unit_mmeasure,
      unit_vmeasure=hitches$funs$unit_vmeasure,
      unit_dmeasure=hitches$funs$unit_dmeasure,
      unit_rmeasure=hitches$funs$unit_rmeasure
    )
  }
)



