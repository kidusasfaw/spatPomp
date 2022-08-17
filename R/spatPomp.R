##' Constructor of the spatPomp object
##'
##' This function constructs a class \sQuote{spatPomp} object, encoding a spatiotemporal partially observed Markov process (\acronym{SpatPOMP}) model together with a uni- or multi-variate time series on a collection of units.
##' Users will typically develop a POMP model for a single unit before embarking on a coupled SpatPOMP analysis.
##' Consequently, we assume some familiarity with \pkg{pomp} and its description by King, Nguyen and Ionides (2016).
##' The \code{spatPomp} class inherits from \code{pomp} with the additional unit structure being a defining feature of the resulting models and inference algorithms.
##'
##' One implements a \acronym{SpatPOMP} model by specifying some or all of its \emph{basic components}, including:
##' \describe{
##' \item{rinit,}{the simulator from the distribution of the latent state process at the zero-time;}
##' \item{rprocess,}{the transition simulator of the latent state process;}
##' \item{dunit_measure,}{the evaluator of the conditional density at a unit's measurement given the unit's latent state;}
##' \item{eunit_measure,}{the evaluator of the expectation of a unit's measurement given the unit's latent state;}
##' \item{munit_measure,}{the evaluator of the moment-matched parameter set given a unit's latent state and some empirical measurement variance;}
##' \item{vunit_measure,}{the evaluator of the variance of a unit's measurement given the unit's latent state;}
##' \item{runit_measure,}{the simulator of a unit's measurement conditional on the unit's latent state;}
##' \item{dprocess,}{the evaluator of the density for transitions of the latent state process;}
##' \item{rmeasure,}{the simulator of the measurements conditional on the latent state;}
##' \item{dmeasure,}{the evaluator of the conditional density of the measurements given the latent state;}
##' \item{rprior,}{the simulator from a prior distribution on the parameters;}
##' \item{dprior,}{the evaluator of the prior density;}
##' \item{skeleton,}{which computes the deterministic skeleton of the unobserved state process;}
##' \item{partrans,}{which performs parameter transformations.}
##' }
##' The basic structure and its rationale are described in Asfaw et al. (2020).
##'
##' Each basic component is supplied via an argument of the same name to \code{spatPomp()}.
##' The five unit-level model components must be provided via C snippets. The remaining components, whose behaviors are inherited from
##' \pkg{pomp} may be furnished using C snippets, \R functions, or pre-compiled native routine available in user-provided dynamically loaded libraries.
##' @param data either a dataframe holding the spatiotemporal data,
##' or an object of class \sQuote{spatPomp}, i.e., the output of another \pkg{spatPomp} calculation.
##' If dataframe, the user must provide the name of the times column using the \code{times} argument and
##' the spatial unit column name using the \code{units} argument. The dataframe provided should be sorted in
##' increasing order of time and unit name respectively, i.e. observation 1 in unit A should come before observation
##' 1 in unit B, which should come before observation 2 in unit A.
##' @param units when \code{data} is a \code{data.frame} this is the name of the column containing the spatial units.
##' @param covar An optional dataframe for supplying covariate information. If provided, there must be two
##' columns that provide the observation time and the observation spatial unit with the same names and arrangement as the \code{data}.
##' @param unit_statenames The names of the components of the latent state. E.g. if the user is constructing an joint SIR model
##' over many spatial units, \code{c('S','I','R')} would be passed.
##' @param unit_accumvars a subset of the \code{unit_statenames} argument that are accumulator variables. See \link[pomp]{accumulator variables}
##' for more on the concept of \pkg{pomp} accumulator variables.
##' @param shared_covarnames If \code{covar} is supplied, covariates that are shared must still be specified for each unit, i.e.,
##' rows with equal values for the same time over all units must be supplied. However, if such covariates exists, supply the names
##' using this argument.
##' @param eunit_measure Evaluator of the expected measurement given the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only C snippets are accepted. The C snippet should assign the scalar approximation to the expected measurement to the pre-defined variable \code{ey} given the latent state and the parameters.
##' For more information, see the examples section below.
##' @param vunit_measure Evaluator of the theoretical measurement variance given the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only C snippets are accepted. The C snippet should assign the scalar approximation to the measurement variance to the pre-defined variable \code{vc} given the latent state and the parameters.
##' For more information, see the examples section below.
##' @param munit_measure Evaluator of a moment-matched parameter set (like the standard deviation parameter of a normal distribution or the size parameter of a negative binomial distribution) given an empirical variance estimate, the latent states and all model parameters.
##' Only Csnippets are accepted. The Csnippet should assign the scalar approximation to the measurement variance parameter to the pre-defined variable corresponding to that parameter, which has been predefined with a \code{M_} prefix. For instance, if the moment-matched parameter is \code{psi}, then the user should assign \code{M_psi} to the moment-matched value.
##' For more information, see the examples section below.
##' @param dunit_measure Evaluator of the unit measurement model density given the measurement, the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar measurement density to the pre-defined variable \code{lik}. The user is encouraged to provide a logged density in an \code{if} condition that checks whether the predefined \code{give_log} variable is true.
##' For more information, see the examples section below.
##' @param runit_measure Simulator of the unit measurement model given the latent states and the model parameters.
##' The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar measurement density to the pre-defined which corresponds to the name of the observation for each unit (e.g. \code{cases} for the measles spatPomp example).
##' For more information, see the examples section below.
##' @param \dots If there are arguments that the user would like to pass to \pkg{pomp}'s basic constructor function's \dots argument,
##' this argument passes them along. Not recommended for this version of \pkg{spatPomp}.
##' @return An object of class \sQuote{spatPomp} representing observations and model components from the spatiotemporal POMP model.
##' @name spatPomp
##' @rdname spatPomp
##'
##' @include pstop.R undefined.R
##' @include rprocess_spec.R safecall.R skeleton_spec.R
##' @include spatPomp_class.R get_covariate_names.R
##' @importFrom rlang .data
##' @inheritParams pomp::pomp
##' @author Kidus Asfaw, Edward L. Ionides, Aaron A. King
##' @references
##'
##' \asfaw2020
##'
##' \king2016
##'
NULL

setGeneric(
  "construct_spatPomp",
  function (data, times, units, ...)
    standardGeneric("construct_spatPomp")
)


##' @rdname spatPomp
##' @export
spatPomp <- function (data, units, times, covar, t0, ...,
  eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure,
  rprocess, rmeasure, dprocess, dmeasure, skeleton, rinit, rprior, dprior,
  unit_statenames, unit_accumvars, shared_covarnames, globals, paramnames, params,
  cdir,cfile, shlib.args, PACKAGE,
  partrans, compile=TRUE, verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("spatPomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (!inherits(data,what=c("data.frame","spatPomp")))
    pStop("spatPomp",sQuote("data")," must be a data frame or an object of ",
      "class ",sQuote("spatPomp"),".")

  ## return as quickly as possible if no work is to be done
  if (is(data,"spatPomp") && missing(times) && missing(t0) &&
        missing(dunit_measure) && missing(eunit_measure) &&
        missing(vunit_measure) && missing(munit_measure) &&
        missing(runit_measure) &&
        missing(rinit) && missing(rprocess) && missing(dprocess) &&
        missing(rmeasure) && missing(dmeasure) && missing(skeleton) &&
        missing(rprior) && missing(dprior) && missing(partrans) &&
        missing(covar) && missing(params) && missing(unit_accumvars) &&
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
      params=params,covar=covar,unit_accumvars=unit_accumvars,
      dunit_measure=dunit_measure,eunit_measure=eunit_measure,
      vunit_measure=vunit_measure,munit_measure=munit_measure,
      runit_measure=runit_measure,unit_statenames=unit_statenames,
      paramnames=paramnames, shared_covarnames=shared_covarnames,PACKAGE=PACKAGE,
      globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args,
      compile=compile, verbose=verbose
    ),
    error = function (e) stop(conditionMessage(e))
  )
}

## Takes care of case where times or units has been set to NULL
setMethod(
  "construct_spatPomp",
  signature=signature(data="ANY", times="ANY", units="ANY"),
  definition = function (data, times, t0, ...) {
    stop(sQuote("times")," should be a single name identifying the column of data that represents",
      " the observation times. ", sQuote("units"), " should be likewise for column that represents",
      " the observation units.")
  }
)

setMethod(
  "construct_spatPomp",
  signature=signature(data="data.frame", times="character", units="character"),
  definition = function (data, times, units, t0, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, unit_accumvars, dunit_measure, eunit_measure,
    vunit_measure, munit_measure, runit_measure, unit_statenames,
    paramnames, shared_covarnames, PACKAGE, globals,
    cdir, cfile, shlib.args, compile, verbose) {

    if (anyDuplicated(names(data)))
      stop("names of data variables must be unique.")

    if (missing(t0)) reqd_arg(NULL,"t0")

    tpos <- match(times,names(data),nomatch=0L)
    upos <- match(units,names(data),nomatch=0L)

    if (length(times) != 1 || tpos == 0L)
      stop(sQuote("times")," does not identify a single column of ",
        sQuote("data")," by name.")
    if (length(units) != 1 || upos == 0L)
      stop(sQuote("units")," does not identify a single column of ",
        sQuote("data")," by name.")

    timename <- times
    unitname <- units

    ## units slot contains unique units. unit_names is an "ordering" of units
    unit_names <- unique(data[[upos]]); U <- length(unit_names)

    ## get observation types
    unit_obsnames <- names(data)[-c(upos,tpos)]
    if(missing(unit_statenames)) unit_statenames <- as.character(NULL)
    if (!missing(unit_accumvars)) pomp_accumvars <- paste0(rep(unit_accumvars,each=U),seq_len(U))
    else {
      pomp_accumvars <- NULL
      unit_accumvars <- as.character(NULL)
    }

    ## if missing workhorses, set to default
    if (missing(rinit)) rinit <- NULL
    if (missing(rprocess) || is.null(rprocess)) {
      rprocess <- rproc_plugin()
    }
    if (missing(dprocess)) dprocess <- NULL
    if (missing(rmeasure)) rmeasure <- NULL
    if (missing(dmeasure)) dmeasure <- NULL
    if (missing(dunit_measure)) dunit_measure <- NULL
    if (missing(eunit_measure)) eunit_measure <- NULL
    if (missing(vunit_measure)) vunit_measure <- NULL
    if (missing(munit_measure)) munit_measure <- NULL
    if (missing(runit_measure)) runit_measure <- NULL
    if (missing(skeleton) || is.null(skeleton)) {
      skeleton <- skel_plugin()
    }
    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL
    if (missing(partrans) || is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(params)) params <- numeric(0)
    if (is.list(params)) params <- unlist(params)
    ## Make data into a dataframe that pomp would expect
    tmp <- seq_along(unit_names)
    names(tmp) <- unit_names
    pomp_data <- data %>% dplyr::mutate(ui = tmp[match(data[,unitname], names(tmp))])
    pomp_data <- pomp_data %>% tidyr::gather(unit_obsnames, key = 'obsname', value = 'val') %>% dplyr::arrange_at(c(timename,'obsname','ui'))
    pomp_data <- pomp_data %>% dplyr::mutate(obsname = paste0(.data$obsname,.data$ui)) %>% dplyr::select(-upos) %>% dplyr::select(-.data$ui)
    pomp_data <- pomp_data %>% tidyr::spread(key = .data$obsname, value = .data$val)
    dat_col_order <- vector(length = U*length(unit_obsnames))
    for(oti in seq_along(unit_obsnames)){
      for(i in seq_len(U)){
        dat_col_order[(oti-1)*U + i] = paste0(unit_obsnames[oti], i)
      }
    }
    pomp_data <- pomp_data[, c(timename, dat_col_order)]
    if(!missing(covar)){
      if(timename %in% names(covar)) tcovar <- timename
      else{
        stop(sQuote("covariate"), ' data.frame should have a time column with the same name as the ',
          'time column of the observation data.frame')
      }
    }
    ## make covariates into a dataframe that pomp would expect
    unit_covarnames <- NULL ## could get overwritten soon
    if(missing(shared_covarnames)) shared_covarnames <- NULL
    if(!missing(covar)){
      upos_cov <- match(unitname, names(covar))
      tpos_cov <- match(tcovar, names(covar))
      cov_col_order <- c()
      if(missing(shared_covarnames)) unit_covarnames <- names(covar)[-c(upos_cov, tpos_cov)]
      else {
        pos_shared_cov <- match(shared_covarnames, names(covar))
        unit_covarnames <- names(covar)[-c(upos_cov, tpos_cov, pos_shared_cov)]
        cov_col_order <- c(cov_col_order, shared_covarnames)
      }
      if(length(unit_covarnames) > 0){
        tmp <- seq_along(unit_names)
        names(tmp) <- unit_names
        pomp_covar <- covar %>% dplyr::mutate(ui = match(covar[,unitname], names(tmp)))
        pomp_covar <- pomp_covar %>% tidyr::gather(unit_covarnames, key = 'covname', value = 'val')
        pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(.data$covname,.data$ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-.data$ui)
        pomp_covar <- pomp_covar %>% tidyr::spread(key = .data$covname, value = .data$val)
        for(cn in unit_covarnames){
          for(i in seq_len(U)){
            cov_col_order = c(cov_col_order, paste0(cn, i))
          }
        }
        pomp_covar <- pomp_covar[, c(timename, cov_col_order)]
        pomp_covar <- pomp::covariate_table(pomp_covar, times=tcovar)
      } else{
        pomp_covar <- pomp::covariate_table(covar, times=tcovar)
      }
    } else {
      pomp_covar <- pomp::covariate_table()
    }

    ## Get all names before call to pomp().
    if(!missing(unit_statenames)) pomp_statenames <- paste0(rep(unit_statenames,each=U),seq_len(U))
    else pomp_statenames <- NULL
    if (!missing(covar)){
      if(missing(shared_covarnames)) pomp_covarnames <- paste0(rep(unit_covarnames,each=U),seq_len(U))
      else {
        if(length(unit_covarnames) == 0) pomp_covarnames <- shared_covarnames
        else pomp_covarnames <- c(shared_covarnames, paste0(rep(unit_covarnames,each=U),seq_len(U)))
      }
    }
    else pomp_covarnames <- NULL
    if (missing(paramnames)) paramnames <- NULL
    if (!missing(paramnames)) mparamnames <- paste("M_", paramnames, sep = "")

    ## We will always have a global giving us the number of spatial units
    if(missing(globals)) globals <- Csnippet(paste0("const int U = ",length(unit_names),";\n"))
    else globals <- Csnippet(paste0(paste0("\nconst int U = ",length(unit_names),";\n"),globals@text))

    ## create the pomp object
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
      accumvars=pomp_accumvars,
      covar = pomp_covar,
      paramnames = paramnames,
      globals = globals,
      cdir = cdir,
      cfile = cfile,
      shlib.args = shlib.args,
      partrans = partrans,
      ...,
      verbose=verbose
    )

    ## Hitch the spatPomp components
    hitches <- pomp::hitch(
                       eunit_measure=eunit_measure,
                       munit_measure=munit_measure,
                       vunit_measure=vunit_measure,
                       dunit_measure=dunit_measure,
                       runit_measure=runit_measure,
                       templates=eval(spatPomp_workhorse_templates),
                       obsnames = paste0(unit_obsnames,"1"),
                       statenames = paste0(unit_statenames,"1"),
                       paramnames=paramnames,
                       mparamnames=mparamnames,
                       covarnames=pomp_covarnames,
                       PACKAGE=PACKAGE,
                       globals=globals,
                       cfile=cfile,
                       cdir=cdir,
                       shlib.args=shlib.args,
                       verbose=verbose
                     )

    pomp::solibs(po) <- hitches$lib
    new("spatPomp",po,
      eunit_measure=hitches$funs$eunit_measure,
      munit_measure=hitches$funs$munit_measure,
      vunit_measure=hitches$funs$vunit_measure,
      dunit_measure=hitches$funs$dunit_measure,
      runit_measure=hitches$funs$runit_measure,
      unit_names=unit_names,
      unit_statenames=unit_statenames,
      unit_accumvars=unit_accumvars,
      unit_obsnames=unit_obsnames,
      unitname=unitname,
      unit_covarnames=as.character(unit_covarnames),
      shared_covarnames=as.character(shared_covarnames))
  }
)

setMethod(
  "construct_spatPomp",
  signature=signature(data="spatPomp", times="NULL", units="NULL"),
  definition = function (data, times, units, t0, timename, unitname, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, paramnames, unit_statenames, covar, shared_covarnames, unit_accumvars,
    dunit_measure, eunit_measure, vunit_measure, munit_measure, runit_measure,
    globals, verbose, PACKAGE, cfile, cdir, shlib.args) {
    times <- data@times
    unit_names <- data@unit_names; U <- length(unit_names)
    if(missing(unit_statenames)) unit_statenames <- data@unit_statenames
    if(length(unit_statenames) == 0) pomp_statenames <- NULL
    else pomp_statenames <- paste0(rep(unit_statenames,each=U),seq_len(U))
    unit_obsnames <- data@unit_obsnames
    if(missing(timename)) timename <- data@timename
    else timename <- as.character(timename)
    if(missing(unitname)) unitname <- data@unitname
    else unitname <- as.character(unitname)

    unit_covarnames <- data@unit_covarnames
    if(missing(shared_covarnames))  shared_covarnames <- data@shared_covarnames
    if(missing(unit_accumvars)) unit_accumvars <- data@unit_accumvars
    if(!missing(covar)){
      if(timename %in% names(covar)) tcovar <- timename
      else{
        stop(sQuote("covariate"), ' data.frame should have a time column with the same name as the ',
          'observation data')
      }
      upos_cov <- match(unitname, names(covar))
      tpos_cov <- match(tcovar, names(covar))
      cov_col_order <- c()
      if(missing(shared_covarnames)) unit_covarnames <- names(covar)[-c(upos_cov, tpos_cov)]
      else {
        pos_shared_cov <- match(shared_covarnames, names(covar))
        unit_covarnames <- names(covar)[-c(upos_cov, tpos_cov, pos_shared_cov)]
        cov_col_order <- c(cov_col_order, shared_covarnames)
      }
      if(length(unit_covarnames) > 0){
        tmp <- seq_along(unit_names)
        names(tmp) <- unit_names
        pomp_covar <- covar %>% dplyr::mutate(ui = match(covar[,unitname], names(tmp)))
        pomp_covar <- pomp_covar %>% tidyr::gather(unit_covarnames, key = 'covname', value = 'val')
        pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(.data$covname,.data$ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-.data$ui)
        pomp_covar <- pomp_covar %>% tidyr::spread(key = .data$covname, value = .data$val)
        for(cn in unit_covarnames){
          for(i in seq_len(U)){
            cov_col_order = c(cov_col_order, paste0(cn, i))
          }
        }
        pomp_covar <- pomp_covar[, c(timename, cov_col_order)]
        pomp_covar <- pomp::covariate_table(pomp_covar, times=tcovar)
      }
    } else pomp_covar <- data@covar

    if (missing(t0)) t0 <- data@t0
    if (missing(rinit)) rinit <- data@rinit
    if (missing(rprocess)) rprocess <- data@rprocess
    else if (is.null(rprocess)) rprocess <- rproc_plugin()
    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure
    if (missing(dunit_measure)) dunit_measure <- data@dunit_measure
    if (missing(munit_measure)) munit_measure <- data@munit_measure
    if (missing(vunit_measure)) vunit_measure <- data@vunit_measure
    if (missing(eunit_measure)) eunit_measure <- data@eunit_measure
    if (missing(runit_measure)) runit_measure <- data@runit_measure
    if (missing(skeleton)) skeleton <- data@skeleton
    else if (is.null(skeleton)) skeleton <- skel_plugin()
    if (missing(rprior)) rprior <- data@rprior
    if (missing(dprior)) dprior <- data@dprior
    if (missing(partrans)) partrans <- data@partrans
    else if (is.null(partrans)) partrans <- parameter_trans()
    if (missing(params) && missing(paramnames)){
      params <- data@params
      paramnames <- names(data@params)
    } else{
      if (!missing(params)) paramnames <- names(params)
    }

    ## Get all names before call to hitch()
    if (!missing(covar)) pomp_covarnames <- paste0(rep(unit_covarnames,each=U),seq_len(U))
    else  pomp_covarnames <- get_covariate_names(data@covar)
    if (!missing(unit_accumvars)) pomp_accumvars <- paste0(rep(unit_accumvars,each=U),seq_len(U))
    else pomp_accumvars <- data@accumvars
    mparamnames <- paste("M_", paramnames, sep = "")

    ## We will always have a global giving us the number of spatial units
    if(missing(globals)) globals <- Csnippet(paste0("const int U = ",length(unit_names),";\n"))
    else globals <- Csnippet(paste0(paste0("\nconst int U = ",length(unit_names),";\n"),globals@text))
    po <- pomp(data = data,
      t0 = t0,
      rprocess = rprocess,
      rmeasure = rmeasure,
      dprocess = dprocess,
      dmeasure = dmeasure,
      skeleton = skeleton,
      rinit = rinit,
      covar = pomp_covar,
      statenames=pomp_statenames,
      accumvars=pomp_accumvars,
      paramnames = paramnames,
      globals = globals,
      cdir = cdir,
      cfile = cfile,
      shlib.args = shlib.args,
      partrans = partrans,
      ...,
      verbose=verbose
    )

    hitches <- pomp::hitch(
                       eunit_measure=eunit_measure,
                       munit_measure=munit_measure,
                       vunit_measure=vunit_measure,
                       dunit_measure=dunit_measure,
                       runit_measure=runit_measure,
                       templates=eval(spatPomp_workhorse_templates),
                       obsnames = paste0(unit_obsnames,"1"),
                       statenames = paste0(unit_statenames,"1"),
                       paramnames=paramnames,
                       mparamnames=mparamnames,
                       covarnames=pomp_covarnames,
                       PACKAGE=PACKAGE,
                       globals=globals,
                       cfile=cfile,
                       cdir=cdir,
                       shlib.args=shlib.args,
                       verbose=verbose
                     )
    pomp::solibs(po) <- hitches$lib
    new(
      "spatPomp",
      po,
      unit_names = unit_names,
      unit_statenames = unit_statenames,
      unit_obsnames = unit_obsnames,
      unitname = unitname,
      shared_covarnames = shared_covarnames,
      unit_covarnames=as.character(unit_covarnames),
      unit_accumvars = unit_accumvars,
      eunit_measure=hitches$funs$eunit_measure,
      munit_measure=hitches$funs$munit_measure,
      vunit_measure=hitches$funs$vunit_measure,
      dunit_measure=hitches$funs$dunit_measure,
      runit_measure=hitches$funs$runit_measure
    )
  }
)
