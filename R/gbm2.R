#' Geometric Brownian motion spatPomp generator
#'
#' Generate a spatPomp object representing a \code{U}-dimensional
#' Geometric Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @return A spatPomp object with the specified dimension and time steps.
#' @examples
#' g <- gbm(U=4, N=20)
#' # See all the model specifications of the object
#' spy(g)
#' @export

gbm <- function(U=5,N=100,delta.t=0.1){   #Create the brownian motion process

U <- U; N <- N; delta.t <- delta.t

#u is first dimension v is second dimension
#u and v go from 1 to U
dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U)) #Distance between 2 units
#Distance related to how correlated the units are (Similar to gravity parameter)
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
gbm_globals <- Csnippet(paste0("#define U ", U, " \n ", dist_C))

#Get observations set up
obs_names <- paste0("Y",1:U)
gbm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)

gbm_unit_statenames <- c("X", "X_bm")
gbm_statenames <- paste0(rep(gbm_unit_statenames, each = U),1:U)

gbm_IVPnames <- paste0(gbm_statenames,"_0") #Initial value of states _0 signifies initial values
gbm_RPnames <- c("rho","sigma","tau")
#Tau= measurement variance of Y | X
  #Fluctuation on what is measured from the process
#Sigma is process noise (standard brownian sigma = 1)
#Brownian = normal (0,t)
gbm_paramnames <- c(gbm_RPnames,gbm_IVPnames)   #Simulate to next time step

gbm_rprocess <- spatPomp_Csnippet("
  double dW[U];
  double pow_rho[U];

//  double X_bm[U]; //Need an X variable to refer to the brownian motion
  //Already passed in as a state name

  int u,v;
//Rho is coupling parameter
//Rho to power of distance is the correlation
  pow_rho[0] = 1;
  for (u=1 ; u < U ; u++) {
    pow_rho[u] = pow_rho[u-1]*rho;
  }

  for (u = 0 ; u < U ; u++) {

    X_bm[u] = log(X[u]);    //Since geometric, the log follows brownian motion

    dW[u] = rnorm(0,sigma*sqrt(dt));
  }
  for (u = 0 ; u < U ; u++) {
    for (v=0; v < U ; v++) {

      X_bm[u] += dW[v]*pow_rho[dist[u][v]]; //Make the brownian motion portion for the log

      //At each time point, own increment as unit and discount other ones surrounding it
      //dw[u]is own noise
      //Things are coupled with new noise term pow-rho
      //X[u] is what already exists
      //When u will iterate for all units
    }
    X[u] = exp(X_bm[u]); //Exponentiate to get geometric form
  }
", unit_statenames = c("X", "X_bm"))

gbm_skel <- spatPomp_Csnippet(" //Process without noise
  //double *X = &X1;
  double *DX = &DX1;
  int u;
  //double dW[U];
  //int u,v;
  for (u = 0 ; u < U ; u++) {
    DX[u] = X[u];
  }
", unit_statenames = c("X"))


gbm_rinit <- Csnippet("
  double *X = &X1;
  const double *X_0 =&X1_0;
  int u;
  for (u = 0; u < U; u++) {
    X[u]=X_0[u];
  }
")

#ASIFIR needs Csnippets that take state
#Evaluates measurement density
gbm_dmeasure <- Csnippet("
  const double *X = &X1;
  const double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  lik=0;
  for (u=0; u<U; u++) lik += dnorm(log(Y[u]),log(X[u]),tau,1)/Y[u];
  if(!give_log) lik = exp(lik) + tol;
")

gbm_unit_emeasure <- Csnippet("
  //Expected based on existing parameters and state that has been simulated y|x
  ey = X* exp(tau*tau/2);
  //X * exp(error)
  //  Turns into X * 1 though
  //https://en.wikipedia.org/wiki/Normal_distribution     MGF formula
")

gbm_unit_mmeasure <- Csnippet("
//Moment matched variance
//Empirical variance: Estimate on process noise based on SSE and SSD
//Bootstrap on variance simulation since don't know true variance
//inverse of vmeasure
//Based on Y|X
//Takes empirical and matches to parameters
  M_tau = sqrt(log(1+sqrt(1+4*(vc/(X*X))))-log(2));
")

gbm_unit_vmeasure <- Csnippet("
//Variance of y | x
//Variance estimated based on parameters
  vc = (X*X) * (exp(2*(tau*tau)) - exp((tau*tau)));
  //https://en.wikipedia.org/wiki/Variance
")

#Two equations of pomp
#Xn from Xn-1 (Brownian Motion)
#Yn from X in (Gaussian N~(0, Tau))

#Brownian motion variance bt is proportional to t
#Increments for brownian motion are independent
#Brownian only continuous process that is Gaussian and stationary

#Draw from measurement distribution
gbm_rmeasure <- Csnippet("
  const double *X = &X1;
  double *Y = &Y1;
  double tol = pow(1.0e-18,U);
  int u;
  for (u=0; u<U; u++) Y[u] = X[u]*exp(rnorm(0,tau+tol)); //Change to make the measure be weighted with process mult bm randomness
")

#Evaluate measurement density
#Compute density Y|X
gbm_unit_dmeasure <- Csnippet("
  double tol = 1.0e-18;
  lik = dnorm(log(Y),log(X),tau,1)/Y +tol;
    //|d/dyg^-1(y)| jacobian adjustment
    //https://www.math.arizona.edu/~jwatkins/f-transform.pdf
  if(!give_log) lik = exp(lik);
")

gbm_unit_rmeasure <- Csnippet("
  double tol = pow(1.0e-18,U);
  double Y;
  Y = X*exp(rnorm(0,tau+tol));
")

#emeasure, mmeasure, and vmeasure needed for ASIFIR and GIRF
gbm_spatPomp <- spatPomp(gbm_data,          #Create the spatPomp Model with the parameters determined from Csnippets
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
               unit_emeasure=gbm_unit_emeasure,
               unit_mmeasure=gbm_unit_mmeasure,
               unit_vmeasure=gbm_unit_vmeasure,
               unit_dmeasure=gbm_unit_dmeasure,
               unit_rmeasure=gbm_unit_rmeasure,
               partrans = parameter_trans(log = c("rho","sigma", "tau")),
               rinit=gbm_rinit
  )


## We need a parameter vector. For now, we initialize the process at zero.
test_ivps <- c(rep(1,U),rep(0,U))  #All initial values 0
                          #2U to accomodate X and X_bm initial value parameters
                          #1 for starting X_i params
                          #0 for starting X_bm params
names(test_ivps) <- gbm_IVPnames
test_params <- c(rho=0.4, sigma=1, tau=1, test_ivps) #Set parameters
      #Tests params, rprocess
simulate(gbm_spatPomp,params=test_params) #Simulate data with the initial parameters
}

