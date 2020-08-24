#' Brownian motion spatPomp generator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional
#' Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @return A spatPomp object with the specified dimension and time steps.
#' @examples
#' b <- bm(U=4, N=20)
#' # See all the model specifications of the object
#' spy(b)
#' @export

ebm <- function(U=5,N=100,rw.sd=0.02,cooling.fraction.50=0.5,delta.t=0.1){

U <- U; N <- N; rw.sd <- rw.sd; cooling.fraction.50 <- cooling.fraction.50; delta.t <- delta.t

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
dist_C <- paste0("const int dist[",U,"][",U,"] = ",dist_C_array,"; ")

# perturbation magnitude
pmag <- vector(length=N)
cooling.fn <- pomp:::mif2.cooling(
  type='geometric',
  fraction=cooling.fraction.50,
  ntimes=N
)
for(i in seq_len(N)){
  pmag[i] <- cooling.fn(i,1)$alpha*rw.sd
}
pmag_C_array <- to_C_array(pmag)
pmag_C <- paste0("const double pmag[",N,"] = ",pmag_C_array,"; ")

ebm_globals <- Csnippet(paste0("#define U ", U, " \n ", dist_C, pmag_C))

obs_names <- paste0("Y",1:U)
ebm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

ebm_unit_statenames <- c("X", "rho", "sigma", "tau")
ebm_statenames <- paste0(rep(ebm_unit_statenames,each=U),1:U)

ebm_IVPnames <- paste0(ebm_statenames,"_0")
ebm_RPnames <- c()
ebm_paramnames <- c(ebm_RPnames,ebm_IVPnames)

ebm_rprocess <- spatPomp_Csnippet("
  double dW[U];
  double pow_rho[U];
  double T_rho, T_sigma, T_tau;
  int u,v;

  T_rho = rnorm(logit(rho[0]),pmag[(int)(ceil(t))]*sqrt(dt));
  T_sigma = rnorm(log(sigma[0]),pmag[(int)(ceil(t))]*sqrt(dt));
  T_tau = rnorm(log(tau[0]),pmag[(int)(ceil(t))]*sqrt(dt));

  for(u = 0; u < U; u++){
    rho[u] = expit(T_rho);
    sigma[u] = exp(T_sigma);
    tau[u] = exp(T_tau);
  }

  pow_rho[0] = 1;
  for (u=1 ; u < U ; u++) {
    pow_rho[u] = pow_rho[u-1]*rho[u];
  }

  for (u = 0 ; u < U ; u++) {
    dW[u] = rnorm(0,sigma[u]*sqrt(dt));
  }
  for (u = 0 ; u < U ; u++) {
    for (v=0; v < U ; v++) {
      X[u] += dW[v]*pow_rho[dist[u][v]];
    }
  }
", unit_statenames = c("X", "rho", "sigma", "tau"))


ebm_rinit <- Csnippet("
  double *X = &X1;
  double *rho = &rho1;
  double *sigma = &sigma1;
  double *tau = &tau1;

  const double *X_0 = &X1_0;
  const double *rho_0 = &rho1_0;
  const double *sigma_0 = &sigma1_0;
  const double *tau_0 = &tau1_0;

  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
    rho[u]=rho_0[u];
    sigma[u]=sigma_0[u];
    tau[u]=tau_0[u];
  }
")


ebm_dmeasure <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  const double *tau = &tau1;

  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau[u],1);
  if(!give_log) lik = exp(lik) + tol;
")

ebm_eunit_measure <- Csnippet("
  ey = X;
")

ebm_vunit_measure <- Csnippet("
  vc = tau*tau;
")

ebm_rmeasure <- Csnippet("
  const double *X = &X1;
  const double *tau = &tau1;

  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau[u]+tol);
")

ebm_dunit_measure <- Csnippet("
  //double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

ebm_runit_measure <- Csnippet("
  double tol = pow(1.0e-18,U);
  double Y;
  Y = rnorm(X,tau+tol);
")

ebm_spatPomp <- spatPomp(ebm_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = ebm_unit_statenames,
               rprocess=euler(ebm_rprocess,delta.t = delta.t),
               paramnames=ebm_paramnames,
               globals=ebm_globals,
               rmeasure=ebm_rmeasure,
               dmeasure=ebm_dmeasure,
               eunit_measure=ebm_eunit_measure,
               vunit_measure=ebm_vunit_measure,
               dunit_measure=ebm_dunit_measure,
               runit_measure=ebm_runit_measure,
               rinit=ebm_rinit
  )


test_ivps <- c(rep(0,U),rep(0.4,U),rep(1,U),rep(1,U))
names(test_ivps) <- ebm_IVPnames
test_params <- c(test_ivps)
simulate(ebm_spatPomp,params=test_params)
}
