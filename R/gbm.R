#' Geometric Brownian motion spatpomp generator
#'
#' Generate a spatpomp object representing a \code{U}-dimensional
#' Geometric Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of time steps to evolve the process.
#' @return A spatpomp object with the specified dimension and time steps.
#' @examples
#' bm(U=4, N=20)
#' @export

gbm <- function(U=5,N=100,delta.t=0.1){

U <- U; N <- N; delta.t <- delta.t

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


obs_names <- paste0("Y",1:U)
gbm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

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
    Xbm[u] = log(X[u])
    dW[u] = rnorm(0,sigma*sqrt(dt));
  }
  for (u = 0 ; u < U ; u++) {
    for (v=0; v < U ; v++) {
      Xbm[u] += dW[v]*pow(rho,dist[u][v]);
    }
    X[u] = exp(Xbm[u])
  }
")

gbm_skel <- Csnippet("
  double *DX = &DX1;
  int u;
  for (u = 0 ; u < U ; u++) {
    DX[u] = X[u];
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
  for (u=0; u<U; u++) lik += dnorm(Y[u],log(X[u]),tau,1);
  if(!give_log) lik = exp(lik) + tol;
")

gbm_rmeasure <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  for (u=0; u<U; u++) Y[u] = rnorm(log(X[u]),tau+tol);
")

gbm_unit_dmeasure <- Csnippet("
  lik = dnorm(Y,log(X),tau,1);
  if(!give_log) lik = exp(lik);
")

gbm_unit_rmeasure <- Csnippet("
  double tol = pow(1.0e-18,U);
  double Y;
  Y = rnorm(log(X),tau+tol);
")

gbm_spatpomp <- spatpomp(gbm_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = gbm_unit_statenames,
               rprocess=euler(gbm_rprocess,delta.t = delta.t),
               skeleton=vectorfield(gbm_skel),
               paramnames=gbm_paramnames,
               globals=gbm_globals,
               rmeasure=gbm_rmeasure,
               dmeasure=gbm_dmeasure,
               unit_dmeasure=gbm_unit_dmeasure,
               unit_rmeasure=gbm_unit_rmeasure,
               partrans = parameter_trans(log = c("rho", "sigma", "tau")),
               rinit=gbm_rinit
  )


## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- rep(0,U)
names(test_ivps) <- gbm_IVPnames
test_params <- c(rho=0.4, sigma=1, tau=2, test_ivps)
simulate(gbm_spatpomp,params=test_params)

}

