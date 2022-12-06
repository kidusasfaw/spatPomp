#' Exact log-likelihood for Brownian motion spatPomp generator
#'
#' Computes the exact likelihood for a model constructed using \code{bm},
#' using the Kalman filter. This model is useful for testing methods
#' in a situation where an exact answer is available
#'
#' @name bm_kalman_logLik
#' @rdname bm_kalman_logLik
#' @author Edward L. Ionides
#' @param bm_object A spatPomp model built using \code{bm}.
#' @param params A parameter vector at which to evaluate the log-likelihood.
#' whereas observations occur every one time unit
#' @return A numeric value for the log-likelihood.
#' @examples
#' # Further examples are provided in the tests directory
#' \dontrun{
#' b <- bm()
#' bm_kalman_logLik(b)
#' }
#' @export
bm_kalman_logLik <- function(bm_object,params=coef(bm_object)){
  U <- length(unit_names(bm_object))
  rows <- matrix(seq(U),nrow=U,ncol=U)
  cols <- matrix(seq(U),nrow=U,ncol=U,byrow=TRUE)
  bm_dmat <- pmin(abs(rows-cols), abs(rows-cols+U), abs(rows-cols-U))
  rootQ <- params["rho"]^bm_dmat * params["sigma"]
  kalmanFilter(bm_object,
    X0=rinit(bm_object),
    A= diag(U),
    Q= rootQ %*% rootQ,
    C=diag(U),
    R=diag(params["tau"]^2, nrow=U)
  )$logLik
}


