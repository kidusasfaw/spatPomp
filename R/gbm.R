#' Geometric Brownian motion spatPomp simulator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional
#' geometric Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time, but
#' an Euler approximation is used for this numerical implementation.
#'
#' @name gbm
#' @rdname gbm
#' @author Kidus Asfaw
#' @family spatPomp model generators
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of time steps to evolve the process.
#' @param delta_t process simulations are performed every \code{delta_t} time units
#' @param delta_obs observations occur every \code{delta_obs} time units
#' @param IVP_values initial value parameters for the latent states
#' @return An object of class \sQuote{spatPomp} representing a simulation from a \code{U}-dimensional
#' geometric Brownian motion
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' g <- gbm(U=4, N=20)
#' # See all the model specifications of the object
#' spy(g)
#' }
#'
#' @references \asfaw2021thesis
#'
#' @export
gbm <- function(U=5, N=100, delta_t=0.1, IVP_values = 1, delta_obs = 1){
  dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
  dmat <- matrix(0,U,U)
  for(u in 1:U) {
    for(v in 1:U) {
      dmat[u,v] <- dist(u,v)
    }
  }
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  dist_C_rows <- apply(dmat,1,to_C_array)
  dist_C_array <- to_C_array(dist_C_rows)
  dist_C <- paste0("const double dist[",U,"][",U,"] = ",dist_C_array,"; ")
  gbm_globals <- Csnippet(paste0("#define U ", U, " \n ", dist_C))

  obs_names <- paste0("Y",(1:U))
  gbm_data <- data.frame(time=rep((1:N)*delta_obs,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

  gbm_unit_statenames <- c("X")
  gbm_statenames <- paste0(gbm_unit_statenames,1:U)

  gbm_IVPnames <- paste0(gbm_statenames,"_0")
  gbm_RPnames <- c("rho","sigma","tau")
  gbm_paramnames <- c(gbm_RPnames,gbm_IVPnames)

  gbm_rprocess <- Csnippet("
    double *X = &X1;
    double Xbm[U];
    double dW[U];
    int u,v;
    for (u = 0 ; u < U ; u++) {
      Xbm[u] = log(X[u]);
      dW[u] = rnorm(0,sigma*sqrt(dt));
    }
    for (u = 0 ; u < U ; u++) {
      for (v=0; v < U ; v++) {
        Xbm[u] += dW[v]*pow(rho,dist[u][v]);
      }
      X[u] = exp(Xbm[u]);
    }
  ")

  gbm_skel <- Csnippet("
  double *DX = &DX1;
  double *X = &X1;
  double cumsigsq[U];
  int u,v;
  for (u = 0; u < U; u++){
     cumsigsq[u] = 0;
  }
  for (u = 0 ; u < U ; u++) {
    for (v = 0 ; v < U ; v++) {
      cumsigsq[u] += pow(sigma*(pow(rho, dist[u][v])), 2);
    }
    DX[u] = X[u]*((cumsigsq[u])/2);
  }
  ")


  gbm_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0=&X1_0;
  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
  }
")


  gbm_dmeasure <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  // Jacobian to get Y|logX since we know logY|logX is Normal(0,tau-squred)
  for (u=0; u<U; u++)lik += (dnorm(log(Y[u]),log(X[u]),tau,1) + log(1/Y[u]));
  if(!give_log) lik = exp(lik) + tol;
  ")

  gbm_rmeasure <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  // Y|X = X*exp(eta) where eta is Normal(0, tau-squared)
  for (u=0; u<U; u++) Y[u] = X[u]*exp(rnorm(0,tau+tol));
  ")

  gbm_dunit_measure <- Csnippet("
  lik = dnorm(log(Y),log(X),tau,1) + log(1/Y);
  if(!give_log) lik = exp(lik);
  ")

  gbm_runit_measure <- Csnippet("
  double tol = pow(1.0e-18,U);
  Y = X*exp(rnorm(0,tau+tol));
  ")

  gbm_eunit_measure <- Csnippet("
  ey = X*exp(tau*tau/2);
  ")

  gbm_munit_measure <- Csnippet("
  M_tau = sqrt(log(0.5 + 0.5*sqrt(1 + (4*vc/(X*X)))));
  ")

  gbm_vunit_measure <- Csnippet("
  vc = X*X*(exp(2*tau*tau) - exp(tau*tau));
  ")

  gbm_spatPomp <- spatPomp(gbm_data,
                           times="time",
                           t0=0,
                           units="unit",
                           unit_statenames = gbm_unit_statenames,
                           rprocess=euler(gbm_rprocess,delta.t = delta_t),
                           skeleton=vectorfield(gbm_skel),
                           paramnames=gbm_paramnames,
                           globals=gbm_globals,
                           rmeasure=gbm_rmeasure,
                           dmeasure=gbm_dmeasure,
                           dunit_measure=gbm_dunit_measure,
                           runit_measure=gbm_runit_measure,
                           eunit_measure=gbm_eunit_measure,
                           munit_measure=gbm_munit_measure,
                           vunit_measure=gbm_vunit_measure,
                           partrans = parameter_trans(logit = c("rho"), log = c("sigma", "tau")),
                           rinit=gbm_rinit
  )
  ## We need a parameter vector. For now, we initialize the Brownian process at zero.
  test_ivps <- rep(IVP_values,U)
  names(test_ivps) <- gbm_IVPnames
  test_params <- c(rho=0.4, sigma=1, tau=1, test_ivps)
  simulate(gbm_spatPomp,params=test_params)

}
