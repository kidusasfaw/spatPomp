#' Calculated log-ARMA log-likelihood benchmark for spatPomp models
#'
#' Fits independent log-ARMA models for each unit, and calculates the conditional
#' log-likelihood for each observation, as well as log-likelihood for
#' each unit and total log-likelihood. A simple tool, but one with
#' practical applicability, as demonstrated by King et al (2008) and
#' Wheeler et al (2023). This function is designed for non-negative 
#' data, and adds 1 to each observation to avoid log(0).
#'
#'
#' @name arma_benchmark
#' @rdname arma_benchmark
#' @author Edward L. Ionides
#' @family utilities
#'
#' @param spo  A spatPomp object
#' @param order A triple (p,d,q) for the ARIMA model fitted to the data. It is
#' intended that d=0
#'
#' @references
#'
#' \king2008
#'
#' \wheeler2023
#'
#'
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' m <- he10(U = 5)
#' arma_benchmark(m)
#' }
#' @export

arma_benchmark <- function(spo,order=c(2,0,1)){
  x <- obs(spo)
  fit <- apply(x,1,function(y,order) arima(log(y + 1), order = order),order=order)
  unit <- sapply(fit,function(z)z$loglik)-apply(x,1,function(y) sum(log(y+1)))
  total <- sum(unit)
  cond <- t(sapply(fit,function(f)unlist(
    dnorm(f$resid,mean=0,sd=sqrt(f$sigma2),log=T)
  ))) - log(x+1)
  list(unit=unit,total=total,cond=cond)
}

 