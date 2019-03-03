#' Linear-Gaussian spatpomp generator
#'
#' Generate a spatpomp object representing a \code{U}-dimensional discrete-time 
#' Ornstein-Uhlenbeck process with \code{N} time steps.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of time steps to evolve the process.
#' @return A spatpomp object with the specified dimension and time steps.
#' @examples
#' ou(5, 100)

ou <- function(U=5,N=100){

U <- 5; N <- 100

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
ou_globals <- Csnippet(paste0("#define U ", U, " \n ", dist_C))


obs_names <- paste0("Y",1:U)
ou_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

ou_unit_statenames <- c("X")
ou_statenames <- paste0(ou_unit_statenames,1:U)

ou_IVPnames <- paste0(ou_statenames,"_0")
ou_RPnames <- c("alpha","rho","sigma","tau")
ou_paramnames <- c(ou_RPnames,ou_IVPnames)

ou_rprocess <- Csnippet("
  double *X = &X1;
  double nextX[U];
  int u,v;

  for (u = 0 ; u < U ; u++) {
    nextX[u] = 0;
    for (v=0; v < U ; v++) {
      nextX[u] += X[v]*pow(rho,dist[u][v]);
    }
    nextX[u] = alpha * nextX[u] + rnorm(0,sigma);
  }
  for (u = 0 ; u < U ; u++) {
    X[u] = nextX[u];
  }
")

ou_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0 =&X1_0;
  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
  }
")


ou_dmeasure <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau,1);
  if(!give_log) lik = exp(lik) + tol;
")

ou_rmeasure <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau+tol);
")

ou_unit_dmeasure <- Csnippet("
  double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

ou <- spatpomp(ou_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = ou_unit_statenames,
               rprocess=discrete_time(ou_rprocess),
               paramnames=ou_paramnames,
               globals=ou_globals,
               rmeasure=ou_rmeasure,
               dmeasure=ou_dmeasure,
               unit_dmeasure=ou_unit_dmeasure,
               partrans = parameter_trans(log = c("rho", "sigma", "tau")),
               rinit=ou_rinit
  )

## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- rep(0,U)
names(test_ivps) <- ou_IVPnames
test_params <- c(alpha=0.4, rho=0.4, sigma=1, tau=1, test_ivps)
simulate(pomp(ou),params=test_params)

}
