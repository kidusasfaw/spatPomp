## package description

#' @name spatPomp-package
#' @title Inference for SpatPOMPs (Spatiotemporal Partially Observed Markov Processes)
#'
#' @description The \pkg{spatPomp} package provides facilities for inference
#' on panel data using spatiotemporal partially-observed Markov process
#' (\acronym{SpatPOMP}) models. To do so, it relies on and extends a number
#' of facilities that the \pkg{pomp} package provides for inference on time
#' series data using partially-observed Markov process (\acronym{POMP}) models.
#'
#' The \pkg{spatPomp} package concerns models consisting of a collection
#' of interacting units. The methods in \pkg{spatPomp} may be applicable
#' whether or not these units correspond to spatial locations.
#'
#' @section Data analysis using \pkg{spatPomp}:
#' The first step in using \pkg{spatPomp} is to encode one's model(s) and data
#'  in objects of class \code{spatPomp}.
#' This can be done via a call to the \link[=spatPomp]{spatPomp} constructor
#' function.
#'
#' @section Extending the \pkg{pomp} platform for developing inference tools:
#' \pkg{spatPomp} extends to panel data the general interface to the
#' components of \acronym{POMP} models provided by \pkg{pomp}. In doing so, it
#' contributes to the goal of the \pkg{pomp} project of facilitating the
#' development of new algorithms in an environment where they can be tested
#' and compared on a growing body of models and datasets.
#'
#' @section Documentation:
#' \pkg{spatPomp} is described by Asfaw et al. (2020)
#'
#' @section License:
#' \pkg{spatPomp} is provided under the \acronym{MIT} License.
#'
#' @references \asfaw2020
#' @author Kidus Asfaw, Joonha Park, Allister Ho, Edward Ionides, Aaron A. King
#'
#' @seealso \link[=pomp-package]{pomp package}
#'
#' @keywords models datasets ts
#'
#' @import methods
"_PACKAGE"

#' @import pomp
#' @useDynLib spatPomp, .registration=TRUE
#' @importFrom stats dnorm runif setNames var
#' @importFrom utils tail
#' @importFrom foreach "%dopar%"
NULL        # replacing NULL by "_PACKAGE" results in roxygen2 adding an
            # \alias{} with the package name, conflicting with functions named
            # after the package
