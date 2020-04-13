

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


  obs_names <- paste0("Y",(1:U))
  gbm_data <- data.frame(time=rep((1:N),U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

  gbm_unit_statenames <- c("X","X_bm")
  gbm_statenames <- paste0(gbm_unit_statenames,1:U)

  gbm_IVPnames <- paste0(gbm_statenames,"_0")
  gbm_RPnames <- c("rho","sigma","tau")
  gbm_paramnames <- c(gbm_RPnames,gbm_IVPnames)


gbm_rprocess <- spatPomp_Csnippet("
  double X_bm[U];
  double dW_bm[U];
  double X[U]
  double pow_rho[U];
  int u,v;

  pow_rho[0] = 1;
  for (u=1 ; u < U ; u++) {
    pow_rho[u] = pow_rho[u-1]*rho;
  }

  for (u = 0 ; u < U ; u++) {
    X_bm[u] = log(X[u]);
    dW_bm[u] = rnorm(0,sigma*sqrt(dt));
  }
  for (u = 0 ; u < U ; u++) {
    for (v=0; v < U ; v++) {
      X_bm[u] += dW_bm[v]*pow_rho[dist[u][v]];
    }
    X[u] = exp(X_bm[u])
  }

", unit_statenames = c("X","X_bm"))

gbm_skel <- Csnippet(" //without noise
  double *X = &X1;
  double *DX = &DX1;
  int u;
  for (u = 0 ; u < U ; u++) {
    DX[u] = X[u];
  }
")

gbm_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0 =&X1_0;
  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
  }
")

gbm_dmeasure <- Csnippet("//measurement density
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  // Jacobian to get Y|logX since we know logY|logX is Normal(0,tau-squred)
  for (u=0; u<U; u++) lik += dnorm(log(Y[u]),log(X[u]),tau,1)+ log(1/Y[u]));
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

gbm_unit_emeasure <- Csnippet(" //x on the B scale, so log(Y) given log(X), E(Y) = X, Var = tau^2 E(Y|log(X)) is not necssarily normal,
                                //find E(e^log(Y)) knwo the dis log(Y) ~ N(log(X), tau^2), do transforamtion, Tau is the measurement variance
  ey = exp(log(X)+0.5*(tau^2));
")

gbm_unit_mmeasure <- Csnippet(" //once answer for vmeasur, sqrt what we got in vmeasure take the inverse of vmeasure
                                //Moment matched variance
                                //inverse of vmeasure
                                // estimate the variance since we don't know the true variance
  M_tau = sqrt(log(1+sqrt(1+4*(vc/(X*X))))-log(2));
")

gbm_unit_vmeasure <- Csnippet(" //variance of y given x log(Y) ~ N(log(x),tau^2). Var(Y|X) similar caculation as emeasure tau is the variance at the BM scale
  vc = exp(2*log(X)+tau^2)*(exp(tau^2)-1);
")

gbm_unit_dmeasure <- Csnippet(" //Compute the density Y|X
  lik = dnorm(log(Y),log(X),tau,1)+log(1/Y) //1 menas the log density not the density
  if(!give_log) lik = exp(lik);
")

gbm_unit_rmeasure <- Csnippet("
  double tol = pow(1.0e-18,U);
  double Y;
  Y = X*exp(rnorm(0,tau+tol)); //the exp of the noise. meansurment variance. measurment model. Y|X.
")

gbm_spatPomp <- spatPomp(gbm_data,
                         times="time",
                         t0=0,
                         units="unit",
                         unit_statenames = gbm_unit_statenames,
                         rprocess=euler(gbm_rprocess,delta.t = delta.t),
                         skeleton=vectorfield(gbm_skel),
                         paramnames=gbm_paramnames,
                         globals=gbm_globals,
                         unit_emeasure=gbm_unit_emeasure,
                         unit_mmeasure=gbm_unit_mmeasure,
                         unit_vmeasure=gbm_unit_vmeasure,
                         unit_dmeasure=gbm_unit_dmeasure,
                         unit_rmeasure=gbm_unit_rmeasure,
                         rmeasure=gbm_rmeasure,
                         dmeasure=gbm_dmeasure,
                         rinit=gbm_rinit,
                         partrans = parameter_trans(log = c("rho","sigma", "tau")),
)

test_ivps <- rep(1,U)
names(test_ivps) <- gbm_IVPnames
test_params <- c(rho=0.1, sigma=0.1, tau=0.1, test_ivps)
simulate(gbm_spatPomp,params=test_params)

}





