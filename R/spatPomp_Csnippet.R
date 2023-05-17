##' C snippets
##'
##' \code{spatPomp_Csnippet()} is used to provide snippets of \proglang{C}
##' code that specify model components. It functions similarly to \code{Csnippet()} from
##' the \pkg{pomp} package; in fact, the output of \code{spatPomp_Csnippet} is an object
##' of class \code{Csnippet}.  It additionally provides some arguments that allow the user
##' to stay focused on model development in the spatiotemporal context  where
##' model size grows.
##'
##' @name spatPomp_Csnippet
##' @rdname spatPomp_Csnippet
##' @include spatPomp_class.R spatPomp.R
##' @author Kidus Asfaw
##'
##' @param code encodes a component of a spatiotemporal POMP model using \proglang{C} code
##' @param method a character string matching the name of the \code{'spatPomp'}
##' argument which the code is designed to specify. This argument is ignored unless
##' needed to correctly specify the Csnippet.
##' @param unit_statenames a subset of the \code{unit_statenames} slot of
##' the \code{spatPomp} object for which we are writing a model. This argument
##' allows the user to get variables that can be indexed conveniently to update
##' states and measurements in a loop. See examples for more details.
##' @param unit_obsnames a subset of the \code{unit_obsnames} slot of
##' the \code{spatPomp} object for which we are writing a model. This argument
##' allows the user to get variables that can be indexed conveniently to update
##' states and measurements in a loop. See examples for more details.
##' @param unit_covarnames if the model has covariate information for each unit,
##' the names of the covariates for each unit can be supplied to this argument.
##' This allows the user to get variables that can be indexed conveniently to
##' use incorporate the covariate information in a loop. See examples for more
##' details.
##' @param unit_ivpnames This argument is particularly useful when specifying the
##' \code{rinit} model component. The \code{paramnames} argument to the
##' \code{spatPomp()} constructor often has names for initial value
##' parameters for the latent states (e.g. \code{S1_0}, \code{S2_0} for the
##' the quantity of susceptibles at unit 1 and unit 2 at the initial time in an
##' SIR model). By supplying \code{unit_ivpnames}, we can get variables
##' that can be easily indexed to reference the initial value parameters (in
##' the previous example, \code{unit_ivpnames=c('S')} we can get a variable
##' named \code{S_0} that we can index as \code{S_0[0]} and \code{S_0[1]} to
##' refer to \code{S1_0} and \code{S2_0}). See examples for more details.
##' @param unit_paramnames This argument is particularly useful when there
##' are non-initial value parameters that are unit-specific.
##' @param unit_vfnames This argument is particularly useful when specifying the
##' \code{skeleton} model component. For all components of the latent state,
##' the user can assume a variable defining the time-derivative is pre-defined (e.g.
##' \code{DS1} and \code{DS2} for the time-derivative of the quantity of the
##' susceptibles at unit 1 and unit 2 in an SIR model). By supplying
##' \code{unit_vfnames}, we can get variables that can be easily indexed to
##' reference these variables (in the previous example,
##' setting \code{unit_vfnames=c('S')} gets us a variable
##' named \code{DS} that we can index as \code{DS[0]} and \code{DS[1]} to
##' refer to \code{DS1} and \code{DS2}). See examples for more details.
##'
##' @return An object of class \sQuote{Csnippet} which represents a model specification in C code.
##' @examples
##' # Set initial states for Brownian motion
##' bm_rinit <- spatPomp_Csnippet(
##'   method = "rinit",
##'   unit_statenames = c("X"),
##'   unit_ivpnames = c("X"),
##'   code = "
##'     for (int u = 0; u < U; u++) {
##'       X[u]=X_0[u];
##'     }
##'   "
##' )
##' # Skeleton for Brownian motion
##' bm_skel <- spatPomp_Csnippet(
##'   method = "skeleton",
##'   unit_statenames = c("X"),
##'   unit_vfnames = c("X"),
##'   code = "
##'       for (int u = 0 ; u < U ; u++) {
##'         DX[u] = 0;
##'       }
##'   "
##')
NULL

setGeneric("spatPomp_Csnippet", function(code,...)standardGeneric("spatPomp_Csnippet"))

##' @name spatPomp_Csnippet-character
##' @rdname spatPomp_Csnippet
##' @aliases spatPomp_Csnippet,character-method
##' @export
setMethod(
  "spatPomp_Csnippet",
  signature=signature(code="character"),
  definition=function(code, method="", unit_statenames, unit_obsnames, unit_covarnames,
    unit_ivpnames, unit_paramnames, unit_vfnames){
    if(missing(unit_statenames) &&
       missing(unit_obsnames) &&
       missing(unit_covarnames) &&
       missing(unit_ivpnames) &&
       missing(unit_paramnames) &&
       missing(unit_vfnames))
      return(pomp::Csnippet(code))
    sn_inits <- on_inits <- cn_inits <- in_inits <- pn_inits <- vn_inits <- character()
      
    if(!missing(unit_statenames)){
      const <- ""
      if(method%in%c("dprocess","rmeasure","dmeasure","skeleton")) const <- "const"
      sn_inits_lhs <- paste0(const, " double *",unit_statenames)
      sn_inits_rhs <- paste0("&", unit_statenames,"1;")
      sn_inits_vec <- paste(sn_inits_lhs, sn_inits_rhs, sep = " = ")
      sn_inits <- paste0(sn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_obsnames)){
      const <- "const"
      if(method%in%c("rmeasure")) const <- ""
      on_inits_lhs <- paste0(const," double *",unit_obsnames)
      on_inits_rhs <- paste0("&", unit_obsnames,"1;")
      on_inits_vec <- paste(on_inits_lhs, on_inits_rhs, sep = " = ")
      on_inits <- paste0(on_inits_vec, collapse = "\n")
    }
    if(!missing(unit_covarnames)){
      cn_inits_lhs <- paste("const double *",unit_covarnames, sep = "")
      cn_inits_rhs <- paste("&", unit_covarnames,"1;",sep="")
      cn_inits_vec <- paste(cn_inits_lhs, cn_inits_rhs, sep = " = ")
      cn_inits <- paste0(cn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_ivpnames)){
      left_ivpnames = paste(unit_ivpnames, "_0", sep = "")
      in_inits_lhs <- paste("const double *",left_ivpnames, sep = "")
      in_inits_rhs <- paste("&", unit_ivpnames,"1_0;",sep="")
      in_inits_vec <- paste(in_inits_lhs, in_inits_rhs, sep = " = ")
      in_inits <- paste0(in_inits_vec, collapse = "\n")
    }
    if(!missing(unit_paramnames)){
      pn_inits_lhs <- paste("const double *",unit_paramnames, sep = "")
      pn_inits_rhs <- paste("&", unit_paramnames,"1;",sep="")
      pn_inits_vec <- paste(pn_inits_lhs, pn_inits_rhs, sep = " = ")
      pn_inits <- paste0(pn_inits_vec, collapse = "\n")
    }
    if(!missing(unit_vfnames)){
      vn_inits_lhs <- paste("double *D",unit_vfnames, sep = "")
      vn_inits_rhs <- paste("&D", unit_vfnames,"1;",sep="")
      vn_inits_vec <- paste(vn_inits_lhs, vn_inits_rhs, sep = " = ")
      vn_inits <- paste0(vn_inits_vec, collapse = "\n")
    }
    all_inits <- paste(sn_inits, on_inits, cn_inits, in_inits, pn_inits, vn_inits, sep = "\n")
    full_csnippet <- paste(all_inits, code, sep = "\n")
    return(pomp::Csnippet(full_csnippet))
  }
)

