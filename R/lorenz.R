#' Lorenz '96 spatPomp simulator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional stochastic Lorenz '96 process with
#' \code{N} measurements made at times \eqn{t_n = n * delta_obs}, simulated using an Euler method
#' with time increment \code{delta_t}.
#'
#' @param U A length-one numeric signifying the number of spatial units for the process.
#' @param N A length-one numeric signifying the number of observations.
#' @param delta_obs A length-one numeric giving the time between observations.
#' @param delta_t A length-one numeric giving the Euler time step for the numerical solution.
#' @param regular_params A named numeric vector containing the values of the \code{F},
#' \code{sigma} and \code{tau} parameters.
#' \code{F=8} is a common value that causes chaotic behavior.
#' @references \lorenz96
#' @return An object of class \sQuote{spatPomp} representing a simulation from a \code{U}-dimensional
#' Lorenz 96 model
#'
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' l <- lorenz(U=5, N=100, delta_t=0.01, delta_obs=1)
#' # See all the model specifications of the object
#' spy(l)
#' }
#' @export
lorenz <- function(U=5,
  N=100,
  delta_t=0.01,
  delta_obs=0.5,
  regular_params=c(F=8, sigma=1, tau=1)){

  if(U<3.5)stop("Please use U >= 4")
  lorenz_globals <- Csnippet(paste0("#define U ", U, "\n"))
  lorenz_unit_statenames <- "X"
  lorenz_obs_names <- paste0("U",1:U)
  lorenz_data <- data.frame(time=rep((1:N)*delta_obs,U),
    unit=rep(lorenz_obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

  ## initial value parameters
  lorenz_IVPnames <- paste0(lorenz_unit_statenames,1:U,"_0")

  ## regular parameters
  lorenz_RPnames <- c("F","sigma","tau")

  ## all parameters
  lorenz_paramnames <- c(lorenz_RPnames,lorenz_IVPnames)

  ## added a condition to prevent numerical instability when the gradient exceeds 1/delta_t
  lorenz_rprocess <- Csnippet("
    double *X = &X1;
    double dXdt[U];
    int u;

    for (u = 2 ; u < U-1 ; u++) {
      dXdt[u] =  (X[u+1]-X[u-2])*X[u-1] - X[u]+F;
    }
    dXdt[0] = (X[1]-X[U-2])*X[U-1] - X[0]+F;
    dXdt[1] = (X[2]-X[U-1])*X[0] - X[1]+F;
    dXdt[U-1] = (X[0]-X[U-3])*X[U-2] - X[U-1]+F;
    for (u = 0 ; u < U ; u++) {
      if(dXdt[u]> 1/dt) dXdt[u] = 1/dt;
      if(dXdt[u]< -1/dt) dXdt[u] = -1/dt;
      X[u] += dXdt[u]*dt + rnorm(0,sigma*sqrt(dt));
    }
  ")

  lorenz_skel <- Csnippet("
    double *X = &X1;
    double dXdt[U];
    double *DX = &DX1;
    int u;
    for (u = 2 ; u < U-1 ; u++) {
      dXdt[u] =  (X[u+1]-X[u-2])*X[u-1] - X[u]+F;
    }
    dXdt[0] = (X[1]-X[U-2])*X[U-1] - X[0]+F;
    dXdt[1] = (X[2]-X[U-1])*X[0] - X[1]+F;
    dXdt[U-1] = (X[0]-X[U-3])*X[U-2] - X[U-1]+F;
    for (u = 0 ; u < U ; u++) {
      DX[u] = dXdt[u];
    }
  ")

  lorenz_rinit <- Csnippet("
    double *X = &X1;
    const double *X_0 =&X1_0;
    int u;
    for (u = 0; u < U; u++) {
      X[u]=X_0[u];
    }
  ")

  lorenz_dmeasure <- Csnippet("
    const double *X = &X1;
    const double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    lik=0;
    for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau,1);
    if(!give_log) lik = exp(lik) + tol;
  ")

  lorenz_rmeasure <- Csnippet("
    const double *X = &X1;
    double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau+tol);
  ")

  lorenz_eunit_measure <- Csnippet("
    ey = X;
  ")

  lorenz_vunit_measure <- Csnippet("
    vc = tau*tau;
  ")

  lorenz_munit_measure <- Csnippet("
    M_tau = sqrt(vc);
  ")

  lorenz_dunit_measure <- Csnippet("
    double tol = 1.0e-18;
    lik = dnorm(Y,X,tau,1);
    if(!give_log) lik = exp(lik);
  ")

  lorenz_runit_measure <- Csnippet("
    double tol = pow(1.0e-18,U);
    double Y;
    Y = rnorm(X,tau+tol);
  ")

  lorenz <- spatPomp(lorenz_data,
    times="time",
    t0=0,
    units="unit",
    unit_statenames = lorenz_unit_statenames,
    rprocess=euler(lorenz_rprocess,delta.t=delta_t),
    skeleton=vectorfield(lorenz_skel),
    paramnames=lorenz_paramnames,
    globals=lorenz_globals,
    rmeasure=lorenz_rmeasure,
    dmeasure=lorenz_dmeasure,
    eunit_measure=lorenz_eunit_measure,
    munit_measure=lorenz_munit_measure,
    vunit_measure=lorenz_vunit_measure,
    dunit_measure=lorenz_dunit_measure,
    runit_measure=lorenz_runit_measure,
    partrans = parameter_trans(log = c("F", "sigma", "tau")),
    rinit=lorenz_rinit)

  ## We need a parameter vector. For now, we initialize the process at zero,
  ## with a small perturbation for state U.
  test_ivps <- c(rep(0,U-1),0.01)
  names(test_ivps) <- lorenz_IVPnames
  test_params <- c(regular_params, test_ivps)
  simulate(lorenz,params=test_params)
}
