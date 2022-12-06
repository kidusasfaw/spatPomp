#' Exact log-likelihood for Brownian motion spatPomp generator
#' with shared or unit-specific parameters
#'
#' Computes the exact likelihood for a model constructed using \code{bm2},
#' using the Kalman filter. This model is useful for testing methods
#' for models with unit-specific parameters, or method such as ibpf
#' which require a unit-specific extension of shared parameters.
#'
#' @name bm2_kalman_logLik
#' @rdname bm2_kalman_logLik
#' @author Edward L. Ionides
#' @param bm2_object A spatPomp model built using \code{bm2}.
#' @param params A parameter vector at which to evaluate the log-likelihood.
#' whereas observations occur every one time unit
#' @return A numeric value for the log-likelihood.
#' @examples
#' # Further examples are provided in the tests directory
#' \dontrun{
#' b <- bm2()
#' bm2_kalman_logLik(b)
#' }
#' @export
bm2_kalman_logLik <- function(bm2_object,params=coef(bm2_object)){
  U <- length(unit_names(bm2_object))
  if(U==1) stop("bm2 designed for U>1")
  ## work out which parameters are unit-specific
  rho_vec <- if("rho2"%in% names(coef(bm2_object)))
    params[paste0("rho",1:U)] else rep(params["rho1"],U)
  sigma_vec <- if("sigma2"%in% names(coef(bm2_object)))
    params[paste0("sigma",1:U)] else rep(params["sigma1"],U)
  tau_vec <- if("tau2"%in% names(coef(bm2_object)))
    params[paste0("tau",1:U)] else rep(params["tau1"],U)  
  U <- length(unit_names(bm2_object))
  rows <- matrix(seq(U),nrow=U,ncol=U)
  cols <- matrix(seq(U),nrow=U,ncol=U,byrow=TRUE)
  bm_dmat <- pmin(abs(rows-cols), abs(rows-cols+U), abs(rows-cols-U))
  rootQ <- matrix(rep(rho_vec,each=U),
    nrow=U,ncol=U,byrow=TRUE)^bm_dmat %*% diag(sigma_vec)
  kalmanFilter(bm2_object,
    X0=rinit(bm2_object),
    A= diag(U),
    Q= rootQ %*% t(rootQ),
    C=diag(U),
    R=diag(tau_vec)^2
  )$logLik
}
