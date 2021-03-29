library(spatPomp)
context("test creation of spatPomp object using measles data")

U = 6
dt = 2/365
fixed_ivps=TRUE
shared_ivps=TRUE
S_0=0.032
E_0=0.00005
I_0=0.00004

birth_lag <- 3*26  # delay until births hit susceptibles, in biweeks

# pre-vaccine biweekly measles reports for the largest 40 UK cities, sorted by size
measlesUK <- spatPomp::measlesUK
measlesUK$city<-as.character(measlesUK$city)
cities <- unique(measlesUK$city)[1:U]
measles_cases <- measlesUK[measlesUK$city %in% cities,c("year","city","cases")]
measles_cases <- measles_cases[measles_cases$year>1949.99,]
measles_covar <- measlesUK[measlesUK$city %in% cities,c("year","city","pop","births")]
u <- split(measles_covar$births,measles_covar$city)
v <- sapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
measles_covar$lag_birthrate <- as.vector(v[,cities])*26
measles_covar$births<- NULL
measles_covarnames <- paste0(rep(c("pop","lag_birthrate"),each=U),1:U)
measles_unit_covarnames <- c("pop","lag_birthrate")

data(city_data_UK)
# Distance between two points on a sphere radius R
# Adapted from geosphere package, which has been cited in the package
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

if(fixed_ivps && shared_ivps){
  ivps_C <- paste0("const double ", c("S_0_shared", "E_0_shared", "I_0_shared"), " = ", c(S_0, I_0, E_0), collapse= ";\n")
  other_ivps_C <- paste0("const double ", c("S1_0", "E1_0", "I1_0"), " = ", c(S_0, I_0, E_0), collapse= ";\n")
  ivps_C <- paste("const int fixed_ivps = 1", "const int shared_ivps = 1", ivps_C,other_ivps_C, ";", sep = ";\n")
}
measles_globals <- Csnippet(
  paste(v_by_g_C, ivps_C, sep = ";\n")
)

measles_unit_statenames <- c('S','E','I','R','C')
measles_RPnames <- c("alpha","iota","psi","R0","gamma","sigma","sigmaSE","cohort","amplitude","mu","rho","g")

if(fixed_ivps && shared_ivps){
  measles_paramnames <- c(measles_RPnames)
} else{
  measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
  measles_IVPnames <- paste0(measles_statenames[1:(4*U)],"_0")
  measles_paramnames <- c(measles_RPnames,measles_IVPnames)
}

measles_rprocess <- Csnippet('
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  //double *W = &W1;
  double powVec[U];
  //double *Acc = &Acc1;
  const double *pop = &pop1;
  const double *lag_birthrate = &lag_birthrate1;
  int obstime = 0;
  int u,v;
  // obstime variable to be used later. See note on if(obstime) conditional
  //if(fabs(((t-floor(t)) / (2.0/52.0)) - (float)(round((t-floor(t)) / (2.0/52.0)))) < 0.001){
       //obstime = 1;
       //Rprintf("t=%f is an observation time\\n",t);
  //}
  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;

  for (u = 0 ; u < U ; u++) {
    // needed for the Ensemble Kalman filter
    // or other methods making real-valued perturbations to the state
    // reulermultinom requires integer-valued double type for states
    S[u] = S[u]>0 ? floor(S[u]) : 0;
    E[u] = E[u]>0 ? floor(E[u]) : 0;
    I[u] = I[u]>0 ? floor(I[u]) : 0;
    R[u] = R[u]>0 ? floor(R[u]) : 0;

    // pre-computing this saves substantial time
    //powVec[u] = pow(I[u]/pop[u],alpha);
    powVec[u] = I[u]/pop[u];
  }

  // These rates could be inside the u loop if some parameters arent shared between units
  rate[1] = mu;			    // natural S death
  rate[2] = sigma;		  // rate of ending of latent stage
  rate[3] = mu;			    // natural E death
  rate[4] = gamma;		  // recovery
  rate[5] = mu;			    // natural I death

  for (u = 0 ; u < U ; u++) {

    // cohort effect
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
      br = cohort*lag_birthrate[u]/dt + (1-cohort)*lag_birthrate[u];
    else
      br = (1.0-cohort)*lag_birthrate[u];

    // expected force of infection
    if(alpha==1.0 && iota==0.0)
      foi = I[u]/pop[u];
    else
      foi = pow( (I[u]+iota)/pop[u],alpha);
    // we follow Park and Ionides (2019) and raise pop to the alpha power
    // He et al (2010) did not do this.

    for (v=0; v < U ; v++) {
      if(v != u)
        foi += g * v_by_g[u][v] * (powVec[v] - powVec[u]) / pop[u];
    }

    // white noise (extrademographic stochasticity)
    dw = rgammawn(sigmaSE,dt);
    rate[0] = beta*foi*dw/dt;  // stochastic force of infection

    // Poisson births
    births = rpois(br*dt);

    // transitions between classes
    reulermultinom(2,S[u],&rate[0],dt,&trans[0]);
    reulermultinom(2,E[u],&rate[2],dt,&trans[2]);
    reulermultinom(2,I[u],&rate[4],dt,&trans[4]);

    S[u] += births   - trans[0] - trans[1];
    E[u] += trans[0] - trans[2] - trans[3];
    I[u] += trans[2] - trans[4] - trans[5];
    R[u] = pop[u] - S[u] - E[u] - I[u];
    //W[u] += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
    C[u] += trans[4];           // true incidence
   }
')

measles_dmeasure <- Csnippet("
  const double *C = &C1;
  const double *cases = &cases1;
  double m,v;
  double tol = 1e-300;
  double mytol = 1e-5;
  int u;

  lik = 0;
  for (u = 0; u < U; u++) {
    m = rho*(C[u]+mytol);
    v = m*(1.0-rho+psi*psi*m);
    // C < 0 can happen in bootstrap methods such as bootgirf
    if (C < 0) {lik += log(tol);} else {
      if (cases[u] > tol) {
        lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-
          pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
      } else {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
      }
    }
  }
  if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
")

measles_rmeasure <- Csnippet("
  const double *C = &C1;
  double *cases = &cases1;
  double m,v;
  double tol = 1.0e-300;
  int u;

  for (u = 0; u < U; u++) {
    m = rho*(C[u]+tol);
    v = m*(1.0-rho+psi*psi*m);
    cases[u] = rnorm(m,sqrt(v)+tol);
    if (cases[u] > 0.0) {
      cases[u] = nearbyint(cases[u]);
    } else {
      cases[u] = 0.0;
    }
  }
")

measles_dunit_measure <- Csnippet('
  double mytol = 1e-5;
  double m = rho*(C+mytol);
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1e-300;
  // C < 0 can happen in bootstrap methods such as bootgirf
  if (C < 0) {lik = 0;} else {
    if (cases > tol) {
      lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-
        pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
    } else {
      lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
    }
  }
  if(give_log) lik = log(lik);
')

measles_eunit_measure <- Csnippet("
ey = rho*C;
")

measles_vunit_measure <- Csnippet("
//consider adding 1 to the variance for the case C = 0
double mytol = 1e-5;
double m;
m = rho*(C+mytol);
vc = m*(1.0-rho+psi*psi*m);
")

measles_munit_measure <- Csnippet("
double binomial_var;
double m;
double mytol = 1e-5;
m = rho*(C+mytol);
binomial_var = rho*(1-rho)*C;
if(vc > binomial_var) {
  M_psi = sqrt(vc - binomial_var)/m;
}
")

measles_rinit <- Csnippet("
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  //double *W = &W1;
  double m;
  const double *pop = &pop1;
  int u;
  if(fixed_ivps && shared_ivps){
    for (u = 0; u < U; u++) {
      m = (float)(pop[u]);
      S[u] = nearbyint(m*S_0_shared);
      I[u] = nearbyint(m*I_0_shared);
      E[u] = nearbyint(m*E_0_shared);
      R[u] = pop[u]-S[u]-E[u]-I[u];
      //W[u] = 0;
      C[u] = 0;
    }
  } else{
    const double *S_0 = &S1_0;
    const double *E_0 = &E1_0;
    const double *I_0 = &I1_0;
    for (u = 0; u < U; u++) {
      m = (float)(pop[u]);
      S[u] = nearbyint(m*S_0[u]);
      I[u] = nearbyint(m*I_0[u]);
      // Use I[u] and the two relvant rates to
      // compute E[u]. (gamma/sigma)*I[u]
      E[u] = nearbyint((gamma/sigma)*(float)(I[u]));
      R[u] = pop[u]-S[u]-E[u]-I[u];
      //W[u] = 0;
      C[u] = 0;
    }
  }
")

measles_skel <- Csnippet('
  double beta, br, seas, foi;
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  //double *W = &W1;
  double *DS = &DS1;
  double *DE = &DE1;
  double *DI = &DI1;
  double *DR = &DR1;
  double *DC = &DC1;
  //double *DW = &DW1;
  double powVec[U];
  const double *pop = &pop1;
  const double *lag_birthrate = &lag_birthrate1;
  int u,v;
  int obstime = 0;
  //if(fabs(((t-floor(t)) / (2.0/52.0)) - (float)(round((t-floor(t)) / (2.0/52.0)))) < 0.001){
       //obstime = 1;
  //}

  // term-time seasonality
   t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;

  // pre-computing this saves substantial time
  for (u = 0 ; u < U ; u++) {
    powVec[u] = pow(I[u]/pop[u],alpha);
  }

  for (u = 0 ; u < U ; u++) {
    // cannot readily put the cohort effect into a vectorfield for the skeleton
    // therefore, we ignore it here.
    // this is okay as long as the skeleton is being used for short-term forecasts
    //    br = lag_birthrate[u];

    // cohort effect, added back in with cohort arriving over a time interval 0.05yr
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*0.05)
      br = cohort*lag_birthrate[u]/0.05 + (1-cohort)*lag_birthrate[u];
    else
      br = (1.0-cohort)*lag_birthrate[u];

    foi = I[u]/pop[u];
    for (v=0; v < U ; v++) {
      if(v != u)
        foi += g * v_by_g[u][v] * (I[v]/pop[v] - I[u]/pop[u]) / pop[u];
    }

    DS[u] = br - (beta*foi + mu)*S[u];
    DE[u] = beta*foi*S[u] - (sigma+mu)*E[u];
    DI[u] = sigma*E[u] - (gamma+mu)*I[u];
    DR[u] = gamma*I[u] - mu*R[u];
    //DW[u] = 0;
    DC[u] = gamma*I[u];
  }
')

if(fixed_ivps && shared_ivps){
  measles_partrans <- parameter_trans(
    log=c("sigma", "gamma", "sigmaSE", "psi", "R0", "g"),
    logit=c("amplitude", "rho")
  )
} else {
  measles_partrans <- parameter_trans(
    log=c("sigma", "gamma", "sigmaSE", "psi", "R0", "g"),
    logit=c("amplitude", "rho",measles_IVPnames)
  )
}

m_full <- spatPomp(measles_cases,
                   units = "city",
                   times = "year",
                   t0 = min(measles_cases$year)-1/26,
                   unit_statenames = measles_unit_statenames,
                   covar = measles_covar,
                   rprocess=euler(measles_rprocess, delta.t=dt),
                   skeleton=vectorfield(measles_skel),
                   unit_accumvars = c("C"),
                   paramnames=measles_paramnames,
                   partrans=measles_partrans,
                   globals=measles_globals,
                   rinit=measles_rinit,
                   dmeasure=measles_dmeasure,
                   eunit_measure=measles_eunit_measure,
                   munit_measure=measles_munit_measure,
                   vunit_measure=measles_vunit_measure,
                   rmeasure=measles_rmeasure,
                   dunit_measure=measles_dunit_measure
)

m_partial <- spatPomp(measles_cases,
                       units = "city",
                       times = "year",
                       t0 = min(measles_cases$year)-1/26,
                       unit_statenames = measles_unit_statenames,
                       covar = measles_covar,
                       rprocess=euler(measles_rprocess, delta.t=dt),
                       unit_accumvars = c("C"),
                       paramnames=measles_paramnames,
                       globals=measles_globals,
                       rinit=measles_rinit,
                       dmeasure=measles_dmeasure,
                       rmeasure=measles_rmeasure,
                       eunit_measure=measles_eunit_measure,
                       dunit_measure=measles_dunit_measure
)

measles_eunit_measure2 <- Csnippet("
                              ey = rho*C;
                              ")

# swap out old eunit_measure for new one
m_partial2 <- spatPomp(m_partial,
                       paramnames=measles_paramnames,
                       eunit_measure=measles_eunit_measure2)

test_that("spatPomp object is created when all and some model components are provided",{
  expect_s4_class(m_full, "spatPomp")
  expect_s4_class(m_partial, "spatPomp")
})

test_that("model component swap-out works",{
  expect_s4_class(m_partial2, "spatPomp")
})



