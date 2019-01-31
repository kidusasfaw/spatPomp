#' Linear-Gaussian spatpomp generator
#'
#' Generate a spatpomp object representing a \code{D}-dimensional Ornstein-Uhlenbeck process with
#' \code{N} time steps.
#'
#' @param D A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of time steps to evolve the process.
#' @return A spatpomp object with the specified dimension and time steps.
#' @examples
#' ou(5, 100)

ou <- function(D=5,N=100){

D <- 3
N <- 10

dist <- function(d,e,n=D) min(abs(d-e),abs(d-e+D),abs(d-e-D))
dmat <- matrix(0,D,D)
for(d in 1:D) {
  for(e in 1:D) {
    dmat[d,e] <- dist(d,e)
  }
}
to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
dist_C_rows <- apply(dmat,1,to_C_array)
dist_C_array <- to_C_array(dist_C_rows)
dist_C <- paste0("const double dist[",D,"][",D,"] = ",dist_C_array,"; ")
ou_globals <- Csnippet(paste0("#define D ", D, " \n ", dist_C))


obs_names <- paste0("Y",1:D)
ou_data <- data.frame(time=rep(1:N,D),unit=rep(obs_names,each=N),Y=rep(NA,D*N),stringsAsFactors=F)

state_names <- paste0("X",1:D)

## initial value parameters
 ivp_names <- paste0(state_names,"_0")

## regular parameters
 rp_names <- c("alpha","rho","sigma","tau")

## all parameters
param_names <- c(rp_names,ivp_names)

ou_rproc <- Csnippet("
  double *X = &X1;
  double nextX[D];
  int d,e;

  for (d = 0 ; d < D ; d++) {
    nextX[d] = 0;
    for (e=0; e < D ; e++) {
      nextX[d] += X[e]*pow(rho,dist[d][e]);
    }
    nextX[d] = alpha * nextX[d] + rnorm(0,sigma);
  }
  for (d = 0 ; d < D ; d++) {
    X[d] = nextX[d];
  }
")

ou_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0 =&X1_0;
  int d;
  for (d = 0; d < D; d++) {
    X[d]=X_0[d];
  }
")


ou_dmeas <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,D);
  int d;
  lik=0;
  for (d=0; d<D; d++) lik += dnorm(Y[d],X[d],tau,1);
  if(!give_log) lik = exp(lik) + tol;
")

ou_rmeas <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,D);
  int d;
  for (d=0; d<D; d++) Y[d] = rnorm(X[d],tau+tol);
")

ou_udmeas <- Csnippet("
  double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

ou <- spatpomp(ou_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = c('X'),
               rprocess=discrete_time(ou_rproc),
               statenames=state_names,
               paramnames=param_names,
               globals=ou_globals,
               rmeasure=ou_rmeas,
               dmeasure=ou_dmeas,
               unit_dmeasure=ou_udmeas,
               partrans = parameter_trans(log = c("rho", "sigma", "tau")),
               rinit=ou_rinit
  )

## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- rep(0,D)
names(test_ivps) <- ivp_names
test_params <- c(alpha=0.4, rho=0.4, sigma=1, tau=1, test_ivps)
simulate(ou,params=test_params)

}
