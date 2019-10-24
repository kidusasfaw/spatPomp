#' Brownian motion spatPomp generator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional
#' Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of time steps to evolve the process.
#' @return A spatPomp object with the specified dimension and time steps.
#' @examples
#' bm(U=4, N=20)
#' @export

bm <- function(U=5,N=100,delta.t=0.1){

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
bm_globals <- Csnippet(paste0("#define U ", U, " \n ", dist_C))


obs_names <- paste0("Y",1:U)
bm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

bm_unit_statenames <- c("X")
bm_statenames <- paste0(bm_unit_statenames,1:U)

bm_IVPnames <- paste0(bm_statenames,"_0")
bm_RPnames <- c("rho","sigma","tau")
bm_paramnames <- c(bm_RPnames,bm_IVPnames)

bm_rprocess <- spatPomp_Csnippet("
  double dW[U];
  int u,v;

  for (u = 0 ; u < U ; u++) {
    dW[u] = rnorm(0,sigma*sqrt(dt));
  }
  for (u = 0 ; u < U ; u++) {
    for (v=0; v < U ; v++) {
      X[u] += dW[v]*pow(rho,dist[u][v]);
    }
  }
", unit_statenames = c("X"))

bm_skel2 <- Csnippet("
  //double *X = &X1;
  double *DX = &DX1;
  int u;
  //double dW[U];
  //int u,v;
  for (u = 0 ; u < U ; u++) {
    DX[u] = 0;
  }
")


bm_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0 =&X1_0;
  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
  }
")


bm_dmeasure <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau,1);
  if(!give_log) lik = exp(lik) + tol;
")

bm_emeasure <- Csnippet("
  ey = X;
")

bm_mmeasure <- Csnippet("
  M_tau = sqrt(vc);
")

bm_vmeasure <- Csnippet("
  vc = tau*tau;
")

bm_rmeasure <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau+tol);
")

bm_unit_dmeasure <- Csnippet("
  //double tol = 1.0e-18;
  lik = dnorm(Y,X,tau,1);
  if(!give_log) lik = exp(lik);
")

bm_unit_rmeasure <- Csnippet("
  double tol = pow(1.0e-18,U);
  double Y;
  Y = rnorm(X,tau+tol);
")

bm_spatPomp <- spatPomp(bm_data,
               times="time",
               t0=0,
               units="unit",
               unit_statenames = bm_unit_statenames,
               rprocess=euler(bm_rprocess,delta.t = delta.t),
               #rprocess=discrete_time(bm_rprocess),
               #skeleton=map(bm_skel, delta.t=delta.t),
               skeleton=vectorfield(bm_skel2),
               paramnames=bm_paramnames,
               globals=bm_globals,
               rmeasure=bm_rmeasure,
               dmeasure=bm_dmeasure,
               emeasure=bm_emeasure,
               mmeasure=bm_mmeasure,
               vmeasure=bm_vmeasure,
               unit_dmeasure=bm_unit_dmeasure,
               unit_rmeasure=bm_unit_rmeasure,
               partrans = parameter_trans(log = c("rho", "sigma", "tau")),
               rinit=bm_rinit
  )


## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- rep(0,U)
names(test_ivps) <- bm_IVPnames
test_params <- c(rho=0.4, sigma=1, tau=2, test_ivps)
simulate(bm_spatPomp,params=test_params)
}

#' @export
girfd_bm <- function(U=5, N = 10, Np = 100, Nguide = 50){
  b <- bm(U = U, N = N)
  # girfd_spatPomp object creation requirements
  bm_Ninter <- length(spat_units(b))
  bm_lookahead <- 1
  bm_Nguide <- Nguide
  bm_Np <- Np
  bm_tol <- 1e-300

  # Output girfd_spatPomp object
  new(
    "girfd_spatPomp",
    b,
    Ninter=bm_Ninter,
    Nguide=bm_Nguide,
    lookahead=bm_lookahead,
    cond.loglik = array(data=numeric(0),dim=c(0,0)),
    Np = as.integer(bm_Np),
    tol= bm_tol,
    loglik=as.double(NA)
  )
}
#' @export
asifird_bm <- function(U=5,
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
  b <- bm(U = U, N = N)
  # asifird_spatPomp object creation requirements
  bm_Np <- Np
  bm_Ninter <- Ninter
  bm_islands <- islands
  bm_nbhd <- nbhd
  bm_tol <- 1e-300

  # Output girfd_spatPomp object
  new(
    "asifird_spatPomp",
    b,
    Np = as.integer(bm_Np),
    Ninter = as.integer(bm_Ninter),
    islands = as.integer(bm_islands),
    nbhd = bm_nbhd,
    tol= bm_tol,
    loglik=as.double(NA)
  )

}


