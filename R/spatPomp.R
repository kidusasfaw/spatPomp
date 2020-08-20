##' Constructor of the basic spatPomp object
##'
##' This function constructs a \sQuote{spatPomp} object, encoding a spatiotemporal partially observed Markov process (\acronym{SpatPOMP}) model together with a uni- or multi-variate time series on a collection of units.
##' Users will typically develop a POMP model for a single unit before embarking on a coupled SpatPOMP analysis.
##' Consequently, we assume some familiarity with \pkg{pomp} and its description by King, Nguyen and Ionides (2016).
##' The \code{spatPomp} class inherits from \code{pomp} with the additional unit structure being a defining feature of the resulting models and inference algorithms.
##'
##' @param eunit_measure Evaluator of the expected measurement given the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar approximation to the expected measurement to the pre-defined variable \code{ey} given the latent state and the parameters.
##' For more information, see the examples section below.
##' @param vunit_measure Evaluator of the theoretical measurement variance given the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar approximation to the measurement variance to the pre-defined variable \code{vc} given the latent state and the parameters.
##' For more information, see the examples section below.
##' @param munit_measure Evaluator of a moment-matched measurement variance parameter (like the standard deviation parameter of a normal distribution or the size parameter of a negative binomial distribution) given an empirical variance estimate, the latent states and all model parameters.
##' Only Csnippets are accepted. The Csnippet should assign the scalar approximation to the measurement variance parameter to the pre-defined variable corresponding to that parameter, which has been predefined with a \code{M_} prefix. For instance, if the moment-matched parameter is \code{psi}, then the user should assign \code{M_psi} to the moment-matched value.
##' For more information, see the examples section below.
##' @param dunit_measure Evaluator of the unit measurement model density given the measurement, the latent states and model parameters. The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar measurement density to the pre-defined variable \code{lik}. The user is encouraged to provide a logged density in an \code{if} condition that checks whether the predefined \code{give_log} variable is true.
##' For more information, see the examples section below.
##' @param runit_measure Simulator of the unit measurement model given the latent states and the model parameters.
##' The \code{unit} variable is pre-defined, which allows the user to specify differing specifications for each unit using \code{if} conditions.
##' Only Csnippets are accepted. The Csnippet should assign the scalar measurement density to the pre-defined which corresponds to the name of the observation for each unit (e.g. \code{cases} for the measles spatPomp example).
##' For more information, see the examples section below.
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
                      eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, unit_statenames, rprocess, rmeasure,
                      dprocess, dmeasure, skeleton, rinit, cdir,cfile, shlib.args, userdata, PACKAGE,
                      globals, statenames, paramnames, unit_obsnames, accumvars, covarnames, shared_covarnames, unit_covarnames,
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

  if (missing(units) && inherits(data,"spatPomp")) unitname <- data@unitname
  else unitname <- units
  if (missing(dprocess)) dprocess <- NULL
  if (missing(rmeasure)) rmeasure <- NULL
  if (missing(dmeasure)) dmeasure <- NULL

  if (missing(skeleton) || is.null(skeleton)) {
    skeleton <- pomp:::skel_plugin()
  }


  if (missing(partrans) || is.null(partrans)) {
    partrans <- parameter_trans()
  }

  if (missing(eunit_measure) && !inherits(data,"spatPomp")) eunit_measure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("eunit_measure")," not specified")
  if (missing(munit_measure) && !inherits(data,"spatPomp")) munit_measure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("munit_measure")," not specified")
  if (missing(vunit_measure) && !inherits(data,"spatPomp")) vunit_measure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("vunit_measure")," not specified")
  if (missing(dunit_measure) && !inherits(data,"spatPomp")) dunit_measure <- function(y,x,t,params,log=FALSE,d,...)
    stop(sQuote("dunit_measure")," not specified")
  if (missing(runit_measure) && !inherits(data,"spatPomp")) runit_measure <- function(x,t,params,log=FALSE,d,...)
    stop(sQuote("runit_measure")," not specified")

  if (missing(unit_statenames)) unit_statenames <- character(0)
  if (missing(unit_covarnames)) unit_covarnames <- character(0)
  if (missing(shared_covarnames)) shared_covarnames <- character(0)


  if (inherits(data, what = "spatPomp")){
    if(!missing(units) && !missing(unit_statenames) && !missing(unit_obsnames))
      stop(ep,sQuote("spatPomp"), "on an existing object can only be used to swap dunit_measure, runit_measure, eunit_measure or munit_measure",call.=FALSE)
    else{
      if(missing(dunit_measure) && missing(runit_measure) && missing(eunit_measure) && missing(munit_measure)){
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
        # inherit from spatPomp (swapping out pomp slots)
        sp <- new("spatPomp",po,
                  unit_covarnames = data@unit_covarnames,
                  shared_covarnames = data@shared_covarnames,
                  runit_measure = data@runit_measure,
                  dunit_measure = data@dunit_measure,
                  eunit_measure = data@eunit_measure,
                  munit_measure = data@munit_measure,
                  unit_names=unit_names(data),
                  unitname=data@unitname,
                  unit_statenames=data@unit_statenames,
                  unit_obsnames = data@unit_obsnames)
        return(sp)
      } else{
        po <- pomp(data)
        covarnames = rownames(data@covar@table)
        unit_statenames = data@unit_statenames
        unit_obsnames = data@unit_obsnames
        paramnames = names(data@params)
        if(!missing(dunit_measure)){
          ## handle dunit_measure C Snippet
          ud_template <- list(
            dunit_measure=list(
              slotname="dunit_measure",
              Cname="__spatPomp_dunit_measure",
              proto=quote(dunit_measure(y,x,t,d,params,log,...)),
              header="\nvoid __spatPomp_dunit_measure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
                  names=unit_obsnames,
                  cref="__y[__obsindex[{%v%}]+unit-1]"
                ),
                lik=list(
                  names="lik",
                  cref="__lik[0]"
                )
              )
            ))
          hitches <- pomp::hitch(
            dunit_measure=dunit_measure,
            templates=ud_template,
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
          # construct new spatpomp object
          pomp:::solibs(po) <- hitches$lib
          # inherit from spatPomp (swapping out spatPomp slots - dunit_measure)
          sp <- new("spatPomp",po,
                    dunit_measure=hitches$funs$dunit_measure,
                    unit_names=data@unit_names,
                    unitname=data@unitname,
                    unit_statenames=data@unit_statenames,
                    unit_obsnames = data@unit_obsnames)
          return(sp)
        } else{
          if(!missing(runit_measure)){
            ur_template <- list(
              runit_measure=list(
                slotname="runit_measure",
                Cname="__spatPomp_runit_measure",
                proto=quote(runit_measure(x,t,d,params,log,...)),
                header="\nvoid __spatPomp_runit_measure (const double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
              runit_measure=runit_measure,
              templates=ur_template,
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
            # construct new spatpomp object
            pomp:::solibs(po) <- hitches$lib
            # inherit from spatPomp (swapping out spatPomp slots - runit_measure)
            sp <- new("spatPomp",po,
                      runit_measure=hitches$funs$runit_measure,
                      unit_names=data@unit_names,
                      unitname=data@unitname,
                      unit_statenames=data@unit_statenames,
                      unit_obsnames = data@unit_obsnames)
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
    unit_obsnames <- names(data)[-c(upos,tpos)]

    # units slot contains unique units. unit_index is an "ordering" of units
    units <- unique(data[[upos]])
    unit_index <- units
    names(unit_index) <- 1:length(units)

    # make data into a dataframe that pomp would expect
    tmp <- names(unit_index)
    names(tmp) <- unit_index
    pomp_data <- data %>% dplyr::mutate(ui = tmp[match(data[,upos_name], names(tmp))])
    pomp_data <- pomp_data %>% tidyr::gather(unit_obsnames, key = 'obsname', value = 'val') %>% dplyr::arrange(pomp_data[,tpos_name], obsname, ui)
    pomp_data <- pomp_data %>% dplyr::mutate(obsname = paste0(obsname,ui)) %>% dplyr::select(-upos) %>% dplyr::select(-ui)
    pomp_data <- pomp_data %>% tidyr::spread(key = obsname, value = val)
    dat_col_order <- vector(length = length(units)*length(unit_obsnames))
    for(ot in unit_obsnames){
      for(i in 1:length(units)){
        dat_col_order[i] = paste0(ot, i)
      }
    }
    pomp_data <- pomp_data[, c(tpos_name, dat_col_order)]
    # make covariates into a dataframe that pomp would expect
    if(!missing(covar) && !missing(tcovar)){
      upos_cov <- match(upos_name, names(covar))
      tpos_cov <- match(tpos_name, names(covar))
      if(!missing(shared_covarnames)){
        spos_cov <- match(shared_covarnames, names(covar))
        covariate_names <- names(covar)[-c(upos_cov, tpos_cov, spos_cov)]
      } else covariate_names <- names(covar)[-c(upos_cov, tpos_cov)]
      tmp <- names(unit_index)
      names(tmp) <- unit_index
      if(length(covariate_names) >0 ){
        pomp_covar <- covar %>% dplyr::mutate(ui = match(covar[,upos_name], names(tmp)))
        pomp_covar <- pomp_covar %>% tidyr::gather(covariate_names, key = 'covname', value = 'val')
        pomp_covar <- pomp_covar %>% dplyr::mutate(covname = paste0(covname,ui)) %>% dplyr::select(-upos_cov) %>% dplyr::select(-ui)
        pomp_covar <- pomp_covar %>% tidyr::spread(key = covname, value = val)
      } else pomp_covar <- covar
      cov_col_order <- c()
      for(cn in covariate_names){
        for(i in 1:length(units)){
          cov_col_order = c(cov_col_order, paste0(cn, i))
        }
      }
      if(!missing(shared_covarnames)) pomp_covar <- pomp_covar[, c(tpos_name, cov_col_order, shared_covarnames)]
      else pomp_covar <- pomp_covar[, c(tpos_name, cov_col_order)]
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
    ## handle dunit_measure C Snippet
    ud_template <- list(
      vunit_measure=list(
        slotname="vunit_measure",
        Cname="__spatPomp_vunit_measure",
        proto=quote(vunit_measure(x,t,d,params,...)),
        header="\nvoid __spatPomp_vunit_measure (double *__vc, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
      munit_measure=list(
        slotname="munit_measure",
        Cname="__spatPomp_munit_measure",
        proto=quote(munit_measure(x,t,d,params,...)),
        header="\nvoid __spatPomp_munit_measure (double *__pm, const double *__x, const double *__p, const double *__vc, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
      eunit_measure=list(
        slotname="eunit_measure",
        Cname="__spatPomp_eunit_measure",
        proto=quote(eunit_measure(y,x,t,d,params,log,...)),
        header="\nvoid __spatPomp_eunit_measure (double *__ey, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
      dunit_measure=list(
        slotname="dunit_measure",
        Cname="__spatPomp_dunit_measure",
        proto=quote(dunit_measure(y,x,t,d,params,log,...)),
        header="\nvoid __spatPomp_dunit_measure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
            names=unit_obsnames,
            cref="__y[__obsindex[{%v%}]+unit-1]"
          ),
          lik=list(
            names="lik",
            cref="__lik[0]"
          )
        )
      ),
      runit_measure=list(
        slotname="runit_measure",
        Cname="__spatPomp_runit_measure",
        proto=quote(runit_measure(x,t,d,params,log,...)),
        header="\nvoid __spatPomp_runit_measure (const double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t, int unit)\n{\n",
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
      eunit_measure=eunit_measure,
      munit_measure=munit_measure,
      vunit_measure=vunit_measure,
      dunit_measure=dunit_measure,
      runit_measure=runit_measure,
      templates=ud_template,
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
    # data.frame -> spatPomp
    new("spatPomp",po,
        eunit_measure=hitches$funs$eunit_measure,
        munit_measure=hitches$funs$munit_measure,
        vunit_measure=hitches$funs$vunit_measure,
        dunit_measure=hitches$funs$dunit_measure,
        runit_measure=hitches$funs$runit_measure,
        unit_names=units,
        unitname=unitname,
        unit_statenames=unit_statenames,
        unit_obsnames = unit_obsnames,
        unit_covarnames=unit_covarnames,
        shared_covarnames=shared_covarnames)

  }
}
