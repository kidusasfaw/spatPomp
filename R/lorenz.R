#' Lorenz 96 spatPomp generator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional stochastic Lorenz 96 process with
#' \code{N} measurements made at times \eqn{t_n= n dt_{obs}}, simulated using an Euler method
#' with time increment $dt$.
#'
#' @param U A length-one numeric signifying the number of spatial units for the process.
#' @param N A length-one numeric signifying the number of observations.
#' @param dt_obs A length-one numeric giving the time between observations.
#' @param dt A length-one numeric giving the Euler time step for the numerical solution.
#'
#' @return A spatPomp object with the specified dimension and time steps.
#'
#' @examples
#' lorenz(U=5, N=100, dt=0.01, dt_obs=1)
#'
#' @export

lorenz <- function(U=5,N=100,dt=0.01,dt_obs=0.5){

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
  int u;

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

lorenz_unit_emeasure <- Csnippet("
  ey = X;
")

lorenz_unit_vmeasure <- Csnippet("
  vc = tau*tau;
")

lorenz_unit_mmeasure <- Csnippet("
  M_tau = sqrt(vc);
")

lorenz_unit_dmeasure <- Csnippet("
  double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

lorenz <- spatPomp(lorenz_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = lorenz_unit_statenames,
               rprocess=euler(lorenz_rprocess,delta.t=dt),
               skeleton=vectorfield(lorenz_skel),
               statenames=lorenz_statenames,
               paramnames=lorenz_paramnames,
               globals=lorenz_globals,
               rmeasure=lorenz_rmeasure,
               dmeasure=lorenz_dmeasure,
               unit_emeasure=lorenz_unit_emeasure,
               unit_mmeasure=lorenz_unit_mmeasure,
               unit_vmeasure=lorenz_unit_vmeasure,
               unit_dmeasure=lorenz_unit_dmeasure,
               partrans = parameter_trans(log = c("F", "sigma", "tau")),
               rinit=lorenz_rinit
  )

## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- c(rep(0,U-1),0.01)
names(test_ivps) <- lorenz_IVPnames
test_params <- c(F=8, sigma=1, tau=1, test_ivps)
simulate(lorenz,params=test_params)
}

#' @export
girfd_lorenz <- function(U=5, N = 10, Np = 100, Nguide = 50, lookahead = 1){
  l <- lorenz(U = U, N = N)
  # girfd_spatPomp object creation requirements
  lorenz_Ninter <- length(spat_units(l))
  lorenz_lookahead <- lookahead
  lorenz_Nguide <- Nguide
  lorenz_Np <- Np
  lorenz_tol <- 1e-300

  # Output girfd_spatPomp object
  new(
    "girfd_spatPomp",
    l,
    Ninter=lorenz_Ninter,
    Nguide=lorenz_Nguide,
    lookahead=lorenz_lookahead,
    cond.loglik = array(data=numeric(0),dim=c(0,0)),
    Np = as.integer(lorenz_Np),
    tol= lorenz_tol,
    loglik=as.double(NA)
  )
}
#' @export
asifird_lorenz <- function(U=5,
                       N = 10,
                       islands = 50,
                       nbhd = function(object, time, unit){
                         nbhd_list = list()
                         if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
                         if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
                         return(nbhd_list)
                       },
                       Np = 10,
                       Ninter = U){
  l <- lorenz(U = U, N = N)
  # asifird_spatPomp object creation requirements
  lorenz_Np <- Np
  lorenz_Ninter <- Ninter
  lorenz_islands <- islands
  lorenz_nbhd <- nbhd
  lorenz_tol <- 1e-300

  # Output girfd_spatPomp object
  new(
    "asifird_spatPomp",
    l,
    Np = as.integer(lorenz_Np),
    Ninter = as.integer(lorenz_Ninter),
    islands = as.integer(lorenz_islands),
    nbhd = lorenz_nbhd,
    tol= lorenz_tol,
    loglik=as.double(NA)
  )

}

