#' Lorenz 96 spatpomp generator
#'
#' Generate a spatpomp object representing a \code{U}-dimensional stochastic Lorenz 96 process with
#' \code{N} measurements made at times \eqn{t_n= n dt_{obs}}, simulated using an Euler method
#' with time increment $dt$. 
#'
#' @param U A length-one numeric signifying the number of spatial units for the process.
#' @param N A length-one numeric signifying the number of observations.
#' @param dt_obs A length-one numeric giving the time between observations.
#' @param dt A length-one numeric giving the Euler time step for the numerical solution.
#' @return A spatpomp object with the specified dimension and time steps.
#' @examples
#' lorenz(U=5, N=100, dt=0.01, dt_obs=1)

lorenz <- function(U=5,N=100,dt=0.01,dt_obs=1){
 
if(U<3.5)stop("Please use U >= 4")

lorenz_globals <- Csnippet(paste0("#define U ", U, "\n"))

lorenz_unit_statenames <- "X"

lorenz_obs_names <- paste0("Y",1:U)
lorenz_data <- data.frame(time=rep((1:N)*dt_obs,U),
  unit=rep(lorenz_obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

lorenz_state_names <- paste0("X",1:U)

## initial value parameters
lorenz_IVPnames <- paste0(lorenz_unit_statenames,1:U,"_0")

## regular parameters
lorenz_RPnames <- c("F","sigma","tau")

## all parameters
lorenz_paramnames <- c(lorenz_RPnames,lorenz_IVPnames)

lorenz_rprocess <- Csnippet("
  double *X = &X1;
  double dXdt[U];
  int u,e;

  for (u = 2 ; u < U-1 ; u++) {
    dXdt[u] =  (X[u+1]-X[u-2])*X[u-1] - X[u]+F;
  }
  dXdt[0] = (X[1]-X[U-2])*X[U-1] - X[0]+F;
  dXdt[1] = (X[2]-X[U-1])*X[0] - X[1]+F;
  dXdt[U-1] = (X[0]-X[U-3])*X[U-2] - X[U-1]+F;
  for (u = 0 ; u < U ; u++) {
    X[u] += dXdt[u]*dt + rnorm(0,sigma*sqrt(dt));
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

lorenz_unit_dmeasure <- Csnippet("
  double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

lorenz <- spatpomp(lorenz_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = lorenz_unit_statenames,
               rprocess=euler(lorenz_rprocess,delta.t=dt),
               statenames=lorenz_statenames, 
               paramnames=lorenz_paramnames,
               globals=lorenz_globals,
               rmeasure=lorenz_rmeasure,
               dmeasure=lorenz_dmeasure,
               unit_dmeasure=lorenz_unit_dmeasure,
               partrans = parameter_trans(log = c("F", "sigma", "tau")),
               rinit=lorenz_rinit
  )

## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- rep(0,U)
names(test_ivps) <- lorenz_IVPnames
test_params <- c(F=8, rho=0.4, sigma=1, tau=1, test_ivps)
simulate(lorenz,params=test_params)

}

