##' Monte Carlo adjusted profile
##'
##' Given a collection of points maximizing the likelihood over a range
##' of fixed values of a focal parameter, this function constructs
##' a profile likelihood confidence interval accommodating both
##' Monte Carlo error in the profile and statistical uncertainty present
##' in the likelihood function.
##'
##' @param lp a vector of profile likelihood evaluations.
##'
##' @param parameter the corresponding values of the focal parameter.
##'
##' @param confidence the required level of the confidence interval.
##'
##' @param lambda the loess parameter used to smooth the profile.
##'
##' @param Ngrid the number of points to evaluate the smoothed profile.
##'
##' @return mcap returns a list including the smoothed profile,
##' a quadratic approximation, and the constructed confidence interval.
##'
##' @author Edward L. Ionides
##'
##' @export
##'
mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000){
  smooth_fit <- loess(lp ~ parameter,span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weights <- rep(0,length(parameter))
  weights[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weights=weights,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit,
                           list(b = parameter_grid, a = -parameter_grid^2)
         )
       ),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared),
       se_mc=sqrt(se_mc_squared),
       se=sqrt(se_total_squared)
  )
}
