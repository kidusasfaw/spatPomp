
#' Book-keeping functions for working with expanded parameters
#'
#' Iterated block particle filters require shared parameters to be
#' expanded into having a value at each unit. expand_params, contract_params
#' and mean_by_unit provide tools for moving between representations.
#' For a unit-specific expansion of a shared parameter, all the values for
#' different units should be the same, and mean_by_unit ensures this
#' by taking an average.
#'
#' These functions assume that expanded parameters have names ending
#' in "1" through "U", where U is the number of units. Contracted parameters, 
#' meaning any parameter that is not expanded, should have a name ending
#' in "1". This numerical suffix convention is useful for writing model-building 
#' code that allows parameters to be either expanded or contracted.
#'
#' @param params Input parameter vector
#' @param expandedParNames character vector of parameters that are, or
#' should be, expanded. These names should have no numerical suffix 1:U.
#' @param U Number of units
#' @rdname param_formats
#' @aliases contract_params, mean_by_unit, expand_params, param_formats
#' @family utilities
#' @export
expand_params <- function(params, expandedParNames,U){
  expanded <- unlist(lapply(expandedParNames,function(par){
    x <- rep(params[paste0(par,'1')],U)
    names(x) <- paste0(par,1:U)
    x
  }))
  unexpandedParNames <- setdiff(names(params),paste0(expandedParNames,'1'))
  unexpanded <- params[unexpandedParNames]
  c(expanded,unexpanded)
}

#' @rdname param_formats
#' @export
contract_params <- function(params, expandedParNames,U,average=FALSE){
if(0){
p_expanded <- c(a1=0,b1=0,b2=1,b3=2,c1=4,c2=4,c3=4)
params <- p_expanded
expandedParNames="c"
U=3
average=F
}
  expanded <- unlist(lapply(expandedParNames,function(par) params[paste0(par,1:U)]))
  unexpanded <- params[setdiff(names(params),names(expanded))]
  contracted <- unlist(lapply(expandedParNames,function(par){
    x <- params[paste0(par,1:U)]
    if(sd(x)>0 & !average) stop ("cannot contract unequal parameters unless average=TRUE")
    x <- mean(x)
    names(x) <- paste0(par,'1')
    x
  }))
  c(unexpanded,contracted)
}

#' @rdname param_formats
#' @export
mean_by_unit <- function(params,expandedParNames,U){
  for(par in expandedParNames){
    params[paste0(par,1:U)] <- mean(params[paste0(par,1:U)])
  }
  params
}

