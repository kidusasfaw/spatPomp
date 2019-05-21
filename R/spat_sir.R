#' SIR model inspired by UK measles: spatpomp generator
#'
#' Generate a spatpomp object for a coupled SIR process 
#' Based on measles in the top-\code{U} most populous cities in England.
#' Model adapted from He et al. (2010) with gravity transport following Park and Ionides (2019).
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatpomp object.
#' @return A spatpomp object.
#' @examples
#' sir1 <- spat_sir()
spat_sir <- function(U=10,Years=20){

cities <- paste0("city",1:U)
year <- (1:(26*Years))/26
sir_empty_data <- data.frame(city=rep(cities,each=Years*26),year=rep(year,U),cases=NA,
   stringsAsFactors=FALSE)


data(city_data_UK)
# Distance between two points on a sphere radius R
# Adapted from geosphere package
distHaversine <- function (p1, p2, r = 6378137)
{
    toRad <- pi/180
    p1 <- p1 * toRad
    p2 <- p2 * toRad
    p = cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2], as.vector(r))
    dLat <- p[, 4] - p[, 2]
    dLon <- p[, 3] - p[, 1]
    a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[, 4]) *
        sin(dLon/2) * sin(dLon/2)
    a <- pmin(a, 1)
    dist <- 2 * atan2(sqrt(a), sqrt(1 - a)) * p[, 5]
    return(as.vector(dist))
}

long_lat <- city_data_UK[1:U,c("lon","lat")]
dmat <- matrix(0,U,U)
for(u1 in 1:U) {
  for(u2 in 1:U) {
    dmat[u1,u2] <- round(distHaversine(long_lat[u1,],long_lat[u2,]) / 1609.344,1)
  }
}

p <- city_data_UK$meanPop[1:U]
v_by_g <- matrix(0,U,U)
dist_mean <- sum(dmat)/(U*(U-1))
p_mean <- mean(p)
for(u1 in 2:U){
  for(u2 in 1:(u1-1)){
    v_by_g[u1,u2] <- (dist_mean*p[u1]*p[u2]) / (dmat[u1,u2] * p_mean^2)
    v_by_g[u2,u1] <- v_by_g[u1,u2]
  }
}
to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
v_by_g_C_rows <- apply(v_by_g,1,to_C_array)
v_by_g_C_array <- to_C_array(v_by_g_C_rows)
v_by_g_C <- Csnippet(paste0("const double v_by_g[",U,"][",U,"] = ",v_by_g_C_array,"; "))
pop <- Csnippet(paste0("const double pop[",U,"] = ",to_C_array(p),"; "))

sir_globals <- Csnippet(
  paste0("const int U = ",U,"; \n ", v_by_g_C,"; \n ", pop)
)

sir_unit_statenames <- c('S','I','R','C')
sir_statenames <- paste0(rep(sir_unit_statenames,each=U),1:U)
sir_params <- c(S_0=0.065, I_0=0.0005, R_0=0.935, beta=0.7*365, amplitude=0.4, gamma=(1/14)*365, sigma=0.01, mu=1/30, rho=0.6, psi=0.1,g=100)
sir_paramnames <- names(sir_params)

sir_rprocess <- Csnippet("
  double seas, foi, dw;
  double rate[4], trans[4];
  double *S = &S1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  int u,v;
  double day;

  // reduced transmission during school summer holiday
  day = (t-floor(t))*365;
  if ( (day<=200)||(day>=250))
      seas = 1.0+amplitude;
    else
      seas = 1.0-amplitude;

  // transition rates
  rate[1] = gamma;		  // recovery
  rate[2] = mu;			  // natural I death
  rate[3] = mu;			  // natural R death

  for (u = 0 ; u < U ; u++) {
    foi = I[u]/pop[u];
    for (v=0; v < U ; v++) {
      if(v != u)
        foi += g * v_by_g[u][v] * (I[v]/pop[v] - I[u]/pop[u]) / pop[u];
    }
    // white noise (extrademographic stochasticity)
    dw = rgammawn(sigma,dt);
    rate[0] = beta*seas*foi*dw/dt;  // stochastic force of infection

    // transitions between classes
    reulermultinom(1,S[u],&rate[0],dt,&trans[0]);
    reulermultinom(2,I[u],&rate[1],dt,&trans[1]);
    reulermultinom(1,R[u],&rate[3],dt,&trans[3]);

    S[u] += trans[2] + trans[3]  - trans[0];
    I[u] += trans[0] - trans[1] - trans[2];
    R[u] += trans[1] - trans[3];
    C[u] += trans[1];       
  }
")

sir_rinit <- Csnippet("
  double *S = &S1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double m;
  int u;
  for (u = 0; u < U; u++) {
    m = pop[u]/(S_0+I_0+R_0);
    S[u] = nearbyint(m*S_0);
    I[u] = nearbyint(m*I_0);
    R[u] = nearbyint(m*R_0);
    C[u] = 0;
  }
")

sir_dmeasure <- Csnippet("
  const double *C = &C1;
  const double *cases = &cases1;
  double m,v;
  double tol = pow(1.0e-18,U);
  int u;

  lik = 0;
  for (u = 0; u < U; u++) {
    m = rho*C[u];
    v = m*(1.0-rho+psi*psi*m);
    if (cases[u] > 0.0) {
      lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
    } else {
      lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
    }
  }
  if(!give_log) lik = exp(lik);
")

sir_rmeasure <- Csnippet("
  const double *C = &C1;
  double *cases = &cases1;
  double m,v;
  double tol = pow(1.0e-18,U);
  int u;

  for (u = 0; u < U; u++) {
    m = rho*C[u];
    v = m*(1.0-rho+psi*psi*m);
    cases[u] = rnorm(m,sqrt(v)+tol);
    if (cases[u] > 0.0) {
      cases[u] = nearbyint(cases[u]);
    } else {
      cases[u] = 0.0;
    }
  }
")

sir_unit_dmeasure <- Csnippet('
                       double m = rho*C;
                       double v = m*(1.0-rho+psi*psi*m);
                       double tol = 1.0e-18;
                       if (cases > 0.0) {
                         lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                       } else {
                           lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                       }
                       ')

sir_spatpomp <- spatpomp(sir_empty_data,
                    units = "city",
                    times = "year",
                    t0 = 0,
                    unit_statenames = sir_unit_statenames,
                    rprocess=euler(sir_rprocess, delta.t=2/365),
                    accumvars = c(paste0("C",1:U)),
                    paramnames=sir_paramnames,
                    globals=sir_globals,
                    rinit=sir_rinit,
                    dmeasure=sir_dmeasure,
                    rmeasure=sir_rmeasure,
                    unit_dmeasure=sir_unit_dmeasure
)

simulate(sir_spatpomp,params=sir_params)

}
