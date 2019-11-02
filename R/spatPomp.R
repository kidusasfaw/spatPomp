##' Constructor of the basic spatPomp object
##'
##' This function constructs a \sQuote{spatPomp} object, encoding a spatiotemporal partially observed Markov process (\acronym{SpatPOMP}) model together with a uni- or multi-variate time series on a collection of units.
##' Users will typically develop a POMP model for a single unit before embarking on a coupled SpatPOMP analysis.
##' Consequently, we assume some familiarity with \pkg{pomp} and its description by King, Nguyen and Ionides (2016).
##' The \code{spatPomp} class inherits from \code{pomp} with the additional unit structure being a defining feature of the resulting models and inference algorithms.
##'
##'##' @param h A user-provided function taking two named arguments: \code{state.vec} (representing the latent state)
##' and \code{param.vec} (representing a parameter vector for the model). It should return a scalar approximation
##' to the expected observed value given a latent state and parameter vector.
##' For more information, see the examples section below.
##' @param theta.to.v A user-provided function taking two named arguments:
##' \code{meas.mean} (representing the observation mean given a latent state - as computed using the \code{h} function above)
##' and \code{param.vec} (representing a parameter vector for the model). It should return a scalar approximation
##' to the variance of the observed value given a latent state and parameter vector.
##' For more information, see the examples section below.
##' @param v.to.theta A user-provided function taking three named arguments:
##' \code{var} (representing an empirical variance), \code{state.vec} (representing a latent state) and \code{param.vec}
##'  (representing a parameter vector for the model). The function should return a parameter vector having observation
##'   noise consistent with variance \code{var} at latent state \code{state.vec} with other parameters given by \code{param.vec}.
##' @name spatPomp
##' @rdname spatPomp
##'
##' @include spatPomp_class.R
##'
##' @param data either a data frame holding the spatiotemporal data,
##' or an object of class \sQuote{spatPomp},
##' i.e., the output of another \pkg{spatPomp} calculation.
##'
##' @inheritParams pomp::pomp
##'
##' @export
spatPomp <- function (data, units, times, covar, tcovar, t0, ...,
  unit_emeasure, unit_mmeasure, unit_vmeasure, unit_dmeasure, unit_rmeasure, unit_statenames, rprocess, rmeasure,
  dprocess, dmeasure, skeleton, rinit, cdir,cfile, shlib.args, userdata, PACKAGE,
  globals, statenames, paramnames, obstypes, accumvars, covarnames,
  partrans, verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("spatPomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (missing(units) && !inherits(data,"spatPomp"))
    stop(ep,sQuote("units")," is a required argument",call.=FALSE)

  if (missing(times) && !inherits(data,"spatPomp"))
    stop(ep,sQuote("times")," is a required argument",call.=FALSE)

  if (missing(rinit)) rinit <- NULL

  if (missing(rprocess) || is.null(rprocess)) {
    rprocess <- pomp:::rproc_plugin()
  }

  if (missing(dprocess)) dprocess <- NULL
  if (missing(rmeasure)) rmeasure <- NULL
  if (missing(dmeasure)) dmeasure <- NULL

  if (missing(skeleton) || is.null(skeleton)) {
    skeleton <- pomp:::skel_plugin()
  }


  if (missing(partrans) || is.null(partrans)) {
    partrans <- parameter_trans()
  }

  if (missing(unit_emeasure) && !inherits(data,"spatPomp")) unit_emeasure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("unit_emeasure")," not specified")
  if (missing(unit_mmeasure) && !inherits(data,"spatPomp")) unit_mmeasure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("unit_mmeasure")," not specified")
  if (missing(unit_vmeasure) && !inherits(data,"spatPomp")) unit_vmeasure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("unit_vmeasure")," not specified")
  if (missing(unit_dmeasure) && !inherits(data,"spatPomp")) unit_dmeasure <- function(y,x,t,params,log=FALSE,d,...)
    stop(sQuote("unit_dmeasure")," not specified")
  if (missing(unit_rmeasure) && !inherits(data,"spatPomp")) unit_rmeasure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("unit_rmeasure")," not specified")

  if (missing(unit_statenames)) unit_statenames <- character(0)

  if (inherits(data, what = "spatPomp")){
    if(!missing(units) && !missing(unit_statenames) && !missing(obstypes))
      stop(ep,sQuote("spatPomp"), "on an existing object can only be used to swap unit_dmeasure and unit_rmeasure",call.=FALSE)
    else{
      if(missing(unit_dmeasure) && missing(unit_rmeasure)){
        # only pomp components are changing
        print("here")
        print(names(data@params))
        po <- pomp(data = data,
                   times=times,
                   t0 = t0,
                   rprocess = rprocess,
                   rmeasure = rmeasure,
                   dprocess = dprocess,
                   dmeasure = dmeasure,
                   skeleton = skeleton,
                   rinit = rinit,
                   statenames=rownames(data@states),
                   accumvars=accumvars,
                   paramnames = paramnames,
                   globals = globals,
                   cdir = cdir,
                   cfile = cfile,
                   shlib.args = shlib.args,
                   partrans = partrans,
                   ...,
                   verbose=verbose
        )
        sp <- new("spatPomp",po,
                  unit_rmeasure = data@unit_rmeasure,
                  unit_dmeasure = data@unit_dmeasure,
                  units=data@units,
                  unit_statenames=data@unit_statenames,
                  obstypes = data@obstypes)
        return(sp)
      } else{
        po <- pomp(data)
        covarnames = rownames(data@covar@table)
        unit_statenames = data@unit_statenames
        obstypes = data@obstypes
        paramnames = names(data@params)
        if(!missing(unit_dmeasure)){
          ## handle unit_dmeasure C Snippet
          ud_template <- list(
            unit_dmeasure=list(
              slotname="unit_dmeasure",
              Cname="__spatPomp_unit_dmeasure",
              proto=quote(unit_dmeasure(y,x,t,d,params,log,...)),
              header="\nvoid __spatPomp_unit_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
                unit_states=list(
                  names=unit_statenames,
                  cref="__x[__stateindex[{%v%}]+unit-1]"
                ),
                obstyp=list(
                  names=obstypes,
                  cref="__y[__obsindex[{%v%}]+unit-1]"
                ),
                lik=list(
                  names="lik",
                  cref="__lik[0]"
                )
              )
            ))
          hitches <- pomp::hitch(
            unit_dmeasure=unit_dmeasure,
            templates=ud_template,
            obsnames = paste0(obstypes,"1"),
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
          # construct new spatpomp object
          pomp:::solibs(po) <- hitches$lib
          sp <- new("spatPomp",po,
              unit_dmeasure=hitches$funs$unit_dmeasure,
              units=data@units,
              unit_statenames=data@unit_statenames,
              obstypes = data@obstypes)
          return(sp)
          } else{
              if(!missing(unit_rmeasure)){
                ur_template <- list(
                  unit_rmeasure=list(
                    slotname="unit_rmeasure",
                    Cname="__spatPomp_unit_rmeasure",
                    proto=quote(unit_rmeasure(x,t,d,params,log,...)),
                    header="\nvoid __spatPomp_unit_rmeasure (const double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
                      unit_states=list(
                        names=unit_statenames,
                        cref="__x[__stateindex[{%v%}]+unit-1]"
                      )
                    )
                  ))
                hitches <- pomp::hitch(
                  unit_rmeasure=unit_rmeasure,
                  templates=ur_template,
                  obsnames = paste0(obstypes,"1"),
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
                # construct new spatpomp object
                pomp:::solibs(po) <- hitches$lib
                sp <- new("spatPomp",po,
                    unit_rmeasure=hitches$funs$unit_rmeasure,
                    units=data@units,
                    unit_statenames=data@unit_statenames,
                    obstypes = data@obstypes)
                return(sp)
              }
        }
      }
    }
  }
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
    unit_index <- units
    names(unit_index) <- 1:length(units)

    # make data into a dataframe that pomp would expect
    tmp <- names(unit_index)
    names(tmp) <- unit_index
    pomp_data <- data %>% dplyr::mutate(ui = tmp[match(data[,upos_name], names(tmp))])
    pomp_data <- pomp_data %>% tidyr::gather(obstypes, key = 'obsname', value = 'val') %>% dplyr::arrange(pomp_data[,tpos_name], obsname, ui)
    pomp_data <- pomp_data %>% dplyr::mutate(obsname = paste0(obsname,ui)) %>% dplyr::select(-upos) %>% dplyr::select(-ui)
    pomp_data <- pomp_data %>% tidyr::spread(key = obsname, value = val)
    dat_col_order <- vector(length = length(units)*length(obstypes))
    for(ot in obstypes){
      for(i in 1:length(units)){
        dat_col_order[i] = paste0(ot, i)
      }
    }
    pomp_data <- pomp_data[, c(tpos_name, dat_col_order)]
    # make covariates into a dataframe that pomp would expect
    if(!missing(covar) && !missing(tcovar)){
      upos_cov <- match(upos_name, names(covar))
      tpos_cov <- match(tpos_name, names(covar))
      covariate_names <- names(covar)[-c(upos_cov, tpos_cov)]
      tmp <- names(unit_index)
      names(tmp) <- unit_index
      pomp_covar <- covar %>% dplyr::mutate(ui = match(covar[,upos_name], names(tmp)))
      pomp_covar <- pomp_covar %>% tidyr::gather(covariate_names, key = 'covname', value = 'val')
      pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(covname,ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-ui)
      pomp_covar <- pomp_covar %>% tidyr::spread(key = covname, value = val)
      cov_col_order <- c()
      for(cn in covariate_names){
        for(i in 1:length(units)){
          cov_col_order = c(cov_col_order, paste0(cn, i))
        }
      }
      pomp_covar <- pomp_covar[, c(tpos_name, cov_col_order)]
      # construct call of covariate_table function
      call_to_covar = list()
      call_to_covar[[1]] <- as.symbol("covariate_table")
      for(col in names(pomp_covar)){
        call_to_covar$col <- pomp_covar[,col]
      }
      call_to_covar$'times'=tpos_name
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

    # get all combinations of unit statenames and units. Concatenate global statenames
    if(!missing(unit_statenames)){
      pomp_statenames <- paste0(rep(unit_statenames,each=length(units)),1:length(units))
    } else pomp_statenames <- character(0)

    # get the observation names of the pomp dataframe.
    obsnames <- names(pomp_data)[-c(tpos)]
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
    ## preliminary error checking
    if (missing(cdir)) cdir <- NULL
    if (missing(cfile)) cfile <- NULL
    if (missing(shlib.args)) shlib.args <- NULL
    if (missing(userdata)) userdata <- list()
    added.userdata <- list(...)
    if (length(added.userdata)>0) {
      message("In ",sQuote("spatPomp"),
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
    if (missing(accumvars)) accumvars <- character(0)
    statenames <- as.character(statenames)
    paramnames <- as.character(paramnames)
    mparamnames <- paste("M_", paramnames, sep = "")
    accumvars <- as.character(accumvars)

    ## check for duplicate names
    if (anyDuplicated(statenames)) {
      stop("all ",sQuote("statenames")," must be unique", call.=FALSE)
    }
    if (anyDuplicated(paramnames)) {
      stop("all ",sQuote("paramnames")," must be unique", call.=FALSE)
    }
    if (anyDuplicated(accumvars)) {
      stop("all ",sQuote("accumvars")," must be unique", call.=FALSE)
    }
    # arrange covariates
    if (missing(covarnames) || length(covarnames)==0)
      if(!is.null(pomp_covar)) covarnames <- as.character(colnames(pomp_covar[-tpos_cov]))
      else covarnames <- NULL
    if (!all(covarnames%in%colnames(pomp_covar))) {
      missing <- covarnames[!(covarnames%in%colnames(covar))]
      stop("covariate(s) ",
        paste(sapply(missing,sQuote),collapse=","),
        " are not among the columns of ",sQuote("covar"),call.=FALSE)
    }
    ## handle unit_dmeasure C Snippet
    ud_template <- list(
      unit_vmeasure=list(
        slotname="unit_vmeasure",
        Cname="__spatPomp_unit_vmeasure",
        proto=quote(unit_vmeasure(x,t,d,params,...)),
        header="\nvoid __spatPomp_unit_vmeasure (double *__vc, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
          unit_states=list(
            names=unit_statenames,
            cref="__x[__stateindex[{%v%}]+unit]"
          ),
          var=list(
            names="vc",
            cref="__vc[0]"
          )
        )
      ),
      unit_mmeasure=list(
        slotname="unit_mmeasure",
        Cname="__spatPomp_unit_mmeasure",
        proto=quote(unit_mmeasure(x,t,d,params,...)),
        header="\nvoid __spatPomp_unit_mmeasure (double *__pm, const double *__x, const double *__p, const double *__vc, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
        footer="\n}\n\n",
        vars=list(
          params=list(
            names=quote(paramnames),
            cref="__p[__parindex[{%v%}]]"
          ),
          mparams=list(
            names=mparamnames,
            cref="__pm[__parindex[{%v%}]]"
          ),
          covars=list(
            names=quote(covarnames),
            cref="__covars[__covindex[{%v%}]]"
          ),
          unit_states=list(
            names=unit_statenames,
            cref="__x[__stateindex[{%v%}]+unit]"
          ),
          var=list(
            names="vc",
            cref="__vc[0]"
          )
        )
      ),
      unit_emeasure=list(
        slotname="unit_emeasure",
        Cname="__spatPomp_unit_emeasure",
        proto=quote(unit_emeasure(y,x,t,d,params,log,...)),
        header="\nvoid __spatPomp_unit_emeasure (double *__ey, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
          unit_states=list(
            names=unit_statenames,
            cref="__x[__stateindex[{%v%}]+unit]"
          ),
          ey=list(
            names="ey",
            cref="__ey[0]"
          )
        )
      ),
      unit_dmeasure=list(
        slotname="unit_dmeasure",
        Cname="__spatPomp_unit_dmeasure",
        proto=quote(unit_dmeasure(y,x,t,d,params,log,...)),
        header="\nvoid __spatPomp_unit_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
          unit_states=list(
            names=unit_statenames,
            cref="__x[__stateindex[{%v%}]+unit-1]"
          ),
          obstyp=list(
            names=obstypes,
            cref="__y[__obsindex[{%v%}]+unit-1]"
          ),
          lik=list(
            names="lik",
            cref="__lik[0]"
          )
        )
      ),
      unit_rmeasure=list(
        slotname="unit_rmeasure",
        Cname="__spatPomp_unit_rmeasure",
        proto=quote(unit_rmeasure(x,t,d,params,log,...)),
        header="\nvoid __spatPomp_unit_rmeasure (const double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
          unit_states=list(
            names=unit_statenames,
            cref="__x[__stateindex[{%v%}]+unit-1]"
          )
        )
      )
    )
    hitches <- pomp::hitch(
      unit_emeasure=unit_emeasure,
      unit_mmeasure=unit_mmeasure,
      unit_vmeasure=unit_vmeasure,
      unit_dmeasure=unit_dmeasure,
      unit_rmeasure=unit_rmeasure,
      templates=ud_template,
      obsnames = paste0(obstypes,"1"),
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
      obstypes = obstypes)

  }
}
