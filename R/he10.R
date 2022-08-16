
#' Measles in UK: spatPomp generator with shared or unit-specific parameters
#'
#' Generate a spatPomp object for measles adding spatial coupling to
#' The model and data from He et al. (2010) with gravity transport as
#' in Park and Ionides (2020). Other transport models may be added in future.
#' The data in the object matches He et al. (2010). The model matches
#' that analysis in the specific case where there is no coupling and all
#' parameters are unit-specific.
#'
#'
#' The code for this spatPomp has duplication with measles(), but in future
#' the two models may diverge. The measles() spatPomp is a simplified
#' situation useful for testing some methods. However, measles() does not
#' permit unit-specific parameters, which he10() allows. Also, 
#' the structure of this spatPomp is compatible with the spatiotemporal 
#' iterated filtering algorithm ibpf(). This requires shared parameters to
#' be represented with a value for each unit, which should be the same for each
#' unit in a valid model instance but may vary between units while optimizing.
#'
#' @param U A length-one numeric signifying the number of cities to be
#' represented in the spatPomp object. Default U=20 gives all the towns
#' studied by He et al., the 10 largest and 10 selected smaller towns.
#' @importFrom utils data read.csv write.table
#' @param dt a numeric (in unit of years) that is used as the Euler time-increment for simulating measles data.
#' @param Tmax Upper time for the window used to construct the object. The lower time is fixed at 1950.0. The default value matches He et al (2010).
#' @param expandedParNames specifies the names of parameters which take unit-specific values. Remaining parameters take a single, shared value for all units.
#' @param basic_params A candidate parameter vector in the basic format, i.e., no unit-specific parameters or unit-related name extensions.
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled measles POMP model.
#' @references
#'
#' \geosphere
#'
#' @note This function goes through a typical workflow of constructing
#' a typical spatPomp object (1-4 below). This allows the user to have a
#' file that replicates the exercise of model building as well as function
#' that creates a typical nonlinear model in epidemiology in case they want
#' to test a new inference methodology. We purposely do not modularize this
#' function because it is not an operational piece of the package and is
#' instead useful as an example.\cr
#' 1. Getting a measurements data.frame with columns for times,
#'    spatial units and measurements.\cr
#' 2. Getting a covariates data.frame with columns for times,
#'    spatial units and covariate data.\cr
#' 3. Constructing model components (latent state initializer,
#'    latent state transition simulator and measurement model). Depending
#'    on the methods used, the user may have to supply a vectorfield to
#'    be integrated that represents the deterministic skeleton of the latent
#'    process.\cr
#' 4. Bringing all the data and model components together to form a
#'    spatPomp object via a call to spatPomp().
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' m <- he10(U = 5)
#' # See all the model specifications of the object
#' spy(m)
#' }
#' @export

# NOTE: Code was written assuming that there are no fixed unit-specific
# parameters.  That could arise if initial values
# are estimated using some other data.

# NOTE: currently assumes measles is loaded from load("data/twentycities.rda").
# In future, this could be included in the package.

he10 <- function(U=6,dt=2/365, Tmax=1964,
  expandedParNames,
  basic_params =c(
    alpha = 1,
    iota = 0,  
    R0 = 30,
    cohort = 0,
    amplitude = 0.5,
    gamma = 52,
    sigma = 52,
    mu = 0.02,
    sigmaSE = 0.15, 
    rho = 0.5,
    psi = 0.15,
    g = 400,
    S_0 = 0.032, 
    E_0 = 0.00005, 
    I_0 = 0.00004
  )
){

## for debugging
## U=6;dt=2/365; Tmax=1964;extendedParNames=c("R0","g","S_0");basic_params =c(alpha = 1,iota = 0,  R0 = 30,cohort = 0,amplitude = 0.5,gamma = 52,sigma = 52,mu = 0.02,sigmaSE = 0.15, rho = 0.5,psi = 0.15,g = 400,S_0 = 0.032, E_0 = 0.00005, I_0 = 0.00004)
  

  if(U>20) stop("U <= 20")
  if(Tmax>1964) stop("Tmax <= 1964")
  birth_lag <- 4 # delay until births hit susceptibles, in years

  # data used for He et al 2010, following their decision
  # to remove 3 data points

  measles_data <- spatPomp::he10measles
  measles_data$date <- as.Date(measles_data$date)
  measles_data$town <- as.character(measles_data$town)

  demog <-  spatPomp::he10demography
  demog$town <- as.character(demog$town)

  # > measles_data[13769+1:5,]
  #            town       date cases
  # 13770 Liverpool 1955-11-04    10
  # 13771 Liverpool 1955-11-11    25
  # 13772 Liverpool 1955-11-18   116
  # 13773 Liverpool 1955-11-25    17
  # 13774 Liverpool 1955-12-02    18

  measles_data[13772,"cases"] <- NA

  # > measles_data[13949+1:5,]
  #            town       date cases
  # 13950 Liverpool 1959-04-17   143
  # 13951 Liverpool 1959-04-24   115
  # 13952 Liverpool 1959-05-01   450
  # 13953 Liverpool 1959-05-08    96
  # 13954 Liverpool 1959-05-15   157

  measles_data[13952,"cases"] <- NA

  # > measles_data[19551+1:5,]
  #             town       date cases
  # 19552 Nottingham 1961-08-18     6
  # 19553 Nottingham 1961-08-25     7
  # 19554 Nottingham 1961-09-01    66
  # 19555 Nottingham 1961-09-08     8
  # 19556 Nottingham 1961-09-15     7

  measles_data[19554,"cases"] <- NA

  mean_pop <- sapply(split(demog,demog$town),function(x) mean(x$pop))
  measles_data <- measles_data[order(mean_pop[measles_data$town],
    -as.numeric(measles_data$date),decreasing=T),]
  towns <-names(sort(mean_pop,decreasing=TRUE))[1:U]
  measles_data$year <- as.integer(format(measles_data$date,"%Y"))
  measles_data %>% 
    dplyr::filter(town%in%towns & year>=1950 & year<Tmax+0.01) %>%
    dplyr::mutate(time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950) %>%
    dplyr::filter(time>1950 & time<Tmax) %>%
    dplyr::select(time,town,cases) -> measles_cases
    
  demog[order(mean_pop[demog$town],
    -as.numeric(demog$year),decreasing=T),] %>%
  dplyr::filter(town%in%towns) -> measles_covar
  colnames(measles_covar)[2] <- "time"

  # note: London starts at 1939, others start at 1940
  u <- split(measles_covar$births,measles_covar$town)
  v <- lapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
  
  
  measles_covar$lag_birthrate <- unlist(v[towns])
  measles_covar$births<- NULL
  measles_covarnames <- paste0(rep(c("pop","lag_birthrate"),each=U),1:U)
  measles_unit_covarnames <- c("pop","lag_birthrate")

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

  long_lat <- spatPomp::he10coordinates[
    match(towns,spatPomp::he10coordinates$town),c("long","lat"),drop=FALSE]
  dmat <- matrix(0,U,U)
  for(u1 in 1:U) {
    for(u2 in 1:U) {
      dmat[u1,u2] <- round(
        distHaversine(long_lat[u1,],long_lat[u2,]) / 1609.344, 1)
    }
  }

  p <- mean_pop[towns]
  v_by_g <- matrix(0,U,U)
  dist_mean <- sum(dmat)/(U*(U-1))
  p_mean <- mean(p)
  if(U>1.5){
    for(u1 in 2:U){
      for(u2 in 1:(u1-1)){
        v_by_g[u1,u2] <- (dist_mean*p[u1]*p[u2]) / (dmat[u1,u2] * p_mean^2)
        v_by_g[u2,u1] <- v_by_g[u1,u2]
      }
    }
  }
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  v_by_g_C_rows <- apply(v_by_g,1,to_C_array)
  v_by_g_C_array <- to_C_array(v_by_g_C_rows)
  v_by_g_C <- Csnippet(paste0("const double v_by_g[",U,"][",U,"] = ",v_by_g_C_array,"; "))

  parNames <- names(basic_params)
  fixedParNames <- setdiff(parNames,expandedParNames)
  set_expanded <- Csnippet(paste0("const int ", expandedParNames,
    "_unit = 1;\n", collapse=" "))
  set_fixed <- Csnippet(paste0("const int ", fixedParNames,
    "_unit = 0;\n", collapse=" "))

  measles_globals <- Csnippet(
    paste(v_by_g_C, set_expanded, set_fixed, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  measles_paramnames <- c(
    if(length(fixedParNames)>0){
      paste0(fixedParNames, "1")
    },
    if(length(expandedParNames)>0){
      paste0(rep(expandedParNames, each=U), 1:U)
    }
  )

  unit_statenames <- c('S','E','I','R','C')
  measles_statenames <- paste0(rep(unit_statenames,each=U),1:U)

  measles_rprocess <- Csnippet('
    const double *amplitude=&amplitude1;
    const double *sigmaSE=&sigmaSE1;
    const double *mu=&mu1;
    const double *g=&g1;
    const double *cohort=&cohort1;
    const double *gamma=&gamma1;
    const double *sigma=&sigma1;
    const double *R0=&R01;
    const double *alpha=&alpha1;
    const double *iota=&iota1;
    double br, beta, seas, foi, dw, births;
    double rate[6], trans[6];
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    double powVec[U];
    const double *pop = &pop1;
    const double *lag_birthrate = &lag_birthrate1;
    int u,v;
    double tol = 1e-18;
    double day = (t-floor(t))*365;

    for (u = 0 ; u < U ; u++) {


      // needed for the Ensemble Kalman filter
      // or other methods making real-valued perturbations to the state
      // reulermultinom requires integer-valued double type for states
      S[u] = S[u]>0 ? floor(S[u]) : 0;
      E[u] = E[u]>0 ? floor(E[u]) : 0;
      I[u] = I[u]>0 ? floor(I[u]) : 0;
      R[u] = R[u]>0 ? floor(R[u]) : 0;

      // pre-computing this saves substantial time
      powVec[u] = pow(I[u]/pop[u],alpha[u*alpha_unit]);
      // powVec[u] = I[u]/pop[u];

    }


    for (u = 0 ; u < U ; u++) {

       if ((day>=7&&day<=100) || (day>=115&&day<=199) || (day>=252&&day<=300) || (day>=308&&day<=356))
        seas = 1.0+amplitude[u*amplitude_unit]*0.2411/0.7589;
      else
        seas = 1.0-amplitude[u*amplitude_unit];

      // transmission rate
      //     beta = R0[u*R0_unit]*(gamma[u*gamma_unit]+mu[u*mu_unit])*seas;

      // for compatilibity with he10
      beta = R0[u*R0_unit]*seas*(1-exp(-(gamma[u*gamma_unit]+mu[u*mu_unit])*dt))/dt;

      rate[1] = mu[u*mu_unit];		// natural S death
      rate[2] = sigma[u*sigma_unit];	// rate of ending of latent stage
      rate[3] = mu[u*mu_unit];		// natural E death
      rate[4] = gamma[u*gamma_unit];	// recovery
      rate[5] = mu[u*mu_unit];		// natural I death

      // cohort effect
      if (fabs(day-251.0+tol) < 0.5*dt*365)
        br = cohort[u*cohort_unit]*lag_birthrate[u]/dt + (1-cohort[u*cohort_unit])*lag_birthrate[u];
      else
        br = (1.0-cohort[u*cohort_unit])*lag_birthrate[u];

      // expected force of infection
      if(alpha[u*alpha_unit]==1.0 && iota[u*iota_unit]==0.0)
        foi = I[u]/pop[u];
      else
        foi = pow( (I[u]+iota[u*iota_unit]),alpha[u*alpha_unit])/pop[u];
      // back to the he10 form where we do not put pop within the power
      //  foi = pow( (I[u]+iota[u*iota_unit])/pop[u],alpha[u*alpha_unit]);
      // following Park and Ionides (2019) we raise pop to the alpha power
      // He et al (2010) did not do this.

      for (v=0; v < U ; v++) {
        if(v != u)
          foi += g[u*g_unit] * v_by_g[u][v] * (powVec[v] - powVec[u]) / pop[u];
      }

      // white noise (extrademographic stochasticity)
      dw = rgammawn(sigmaSE[u*sigmaSE_unit],dt);
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
      C[u] += trans[4];           // true incidence
     }
  ')

  measles_dmeasure <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;
    const double *rho=&rho1;
    const double *psi=&psi1;
    double m,v;
    double tol = 1e-300;
    double mytol = 1e-5;

// smaller tolerances can be appropriate for higher dimensions
// to check compatilibity with he10 for single towns:
//    double tol = 1e-18;
//    double mytol = 0;

    int u;

    lik = 0;
    for (u = 0; u < U; u++) {
      m = rho[u*rho_unit]*(C[u]+mytol);
      v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);

      // Deal with NA measurements by omitting them
      if(!(ISNA(cases[u]))){
        // C < 0 can happen in bootstrap methods such as bootgirf
        if (C[u] < 0) {lik += log(tol);} else {
          if (cases[u] > tol) {
            lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-
              pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
          } else {
              lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
          }
        }
      }
    }

    if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
    
  ")

  measles_rmeasure <- Csnippet("
    const double *C = &C1;
    double *cases = &cases1;
    const double *rho = &rho1;
    const double *psi = &psi1;
    double m,v;
    double tol = 1.0e-300;
    int u;

    for (u = 0; u < U; u++) {
      m = rho[u*rho_unit]*(C[u]+tol);
      v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
      cases[u] = rnorm(m,sqrt(v)+tol);
      if (cases[u] > 0.0) {
        cases[u] = nearbyint(cases[u]);
      } else {
        cases[u] = 0.0;
      }
    }
  ")

  measles_dunit_measure <- Csnippet('
    const double *rho = &rho1;
    const double *psi = &psi1;
    double mytol = 1e-5;
    double m = rho[(u-1)*rho_unit]*(C+mytol);
    double v = m*(1.0-rho[(u-1)*rho_unit]+psi[(u-1)*psi_unit]*psi[(u-1)*psi_unit]*m);
    double tol = 1e-300;
    // C < 0 can happen in bootstrap methods such as bootgirf
    if(ISNA(cases)) {lik=1;} else { 
      if (C < 0) {lik = 0;} else {
        if (cases > tol) {
          lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-
            pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
        } else {
          lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
        }
      }
    }
    if(give_log) lik = log(lik);
  ')

  measles_eunit_measure <- Csnippet("
    const double *rho = &rho1;
    ey = rho[u*rho_unit]*C;
  ")

  measles_vunit_measure <- Csnippet("
    //consider adding 1 to the variance for the case C = 0
    const double *rho = &rho1;
    const double *psi = &psi1;
    double mytol = 1e-5;
    double m;
    m = rho[u*rho_unit]*(C+mytol);
    vc = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
  ")

  measles_munit_measure <- Csnippet("
    const double *rho = &rho1;
    double binomial_var;
    double m;
    double mytol = 1e-5;
    m = rho[u*rho_unit]*(C+mytol);
    binomial_var = rho[u*rho_unit]*(1-rho[u*rho_unit])*C;
    if(vc > binomial_var) {
      M_psi1 = sqrt(vc - binomial_var)/m;
    }
  ")

  measles_rinit <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    const double *S_0 = &S_01;
    const double *E_0 = &E_01;
    const double *I_0 = &I_01;
    double m;
    const double *pop = &pop1;
    int u;
    for(u = 0; u < U; u++) {
      m = (float)(pop[u]);
      S[u] = nearbyint(m*S_0[u*S_0_unit]);
      I[u] = nearbyint(m*I_0[u*I_0_unit]);
      E[u] = nearbyint(m*E_0[u*E_0_unit]);
      R[u] = pop[u]-S[u]-E[u]-I[u];
      C[u] = 0;
      // in any practical model fit, we expect R>0
      // though the model does not strictly enforce that 
    }
  ")

  measles_skel <- Csnippet('
    double beta, br, seas, foi;
    const double *amplitude=&amplitude1;
    const double *mu=&mu1;
    const double *g=&g1;
    const double *cohort=&cohort1;
    const double *gamma=&gamma1;
    const double *sigma=&sigma1;
    const double *R0=&R01;
    const double *alpha=&alpha1;
    const double *iota=&iota1;
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    double *DS = &DS1;
    double *DE = &DE1;
    double *DI = &DI1;
    double *DR = &DR1;
    double *DC = &DC1;
    double powVec[U];
    const double *pop = &pop1;
    const double *lag_birthrate = &lag_birthrate1;
    int u,v;

    // pre-computing this saves substantial time
    for (u = 0 ; u < U ; u++) {
      powVec[u] = pow(I[u]/pop[u],alpha[u*alpha_unit]);
    }

    for (u = 0 ; u < U ; u++) {

      // term-time seasonality
      t = (t-floor(t))*365.25;
      if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
        seas = 1.0+amplitude[u*amplitude_unit]*0.2411/0.7589;
      else
        seas = 1.0-amplitude[u*amplitude_unit];

    // transmission rate
    beta = R0[u*R0_unit]*(gamma[u*gamma_unit]+mu[u*amplitude_unit])*seas;

      // cannot readily put the cohort effect into a vectorfield for the skeleton
      // therefore, we ignore it here.
      // this is okay as long as the skeleton is being used for short-term forecasts
      //    br = lag_birthrate[u];

      // cohort effect, added back in with cohort arriving over a time interval 0.05yr
      if (fabs(t-floor(t)-251.0/365.0) < 0.5*0.05)
        br = cohort[u*cohort_unit]*lag_birthrate[u]/0.05 + (1-cohort[u*cohort_unit])*lag_birthrate[u];
      else
        br = (1.0-cohort[u*cohort_unit])*lag_birthrate[u];

      foi = I[u]/pop[u];
      for (v=0; v < U ; v++) {
        if(v != u)
          foi += g[u*g_unit] * v_by_g[u][v] * (I[v]/pop[v] - I[u]/pop[u]) / pop[u];
      }

      DS[u] = br - (beta*foi + mu[u*mu_unit])*S[u];
      DE[u] = beta*foi*S[u] - (sigma[u*sigma_unit]+mu[u*mu_unit])*E[u];
      DI[u] = sigma[u*sigma_unit]*E[u] - (gamma[u*gamma_unit]+mu[u*mu_unit])*I[u];
      DR[u] = gamma[u*gamma_unit]*I[u] - mu[u*mu_unit]*R[u];
      DC[u] = gamma[u*gamma_unit]*I[u];
    }
  ')

basic_log_names <- c("sigma", "gamma", "sigmaSE", "psi", "R0", "g","iota","mu")
basic_log_names <- setdiff(basic_log_names,fixedParNames)

basic_logit_names <- c("amplitude", "alpha","cohort","rho","S_0", "E_0", "I_0")
basic_logit_names <- setdiff(basic_logit_names,fixedParNames)

log_names <- unlist(lapply(basic_log_names, function(x,U) paste0(x,1:U),U))
logit_names <- unlist(lapply(basic_logit_names, function(x,U) paste0(x,1:U),U))
  
## it is possible for S+E+I to be greater than P,
## in which case R is negative, but that is not necessarily a critical problem.

measles_partrans <- parameter_trans(log=log_names,logit=logit_names)

m1 <-  spatPomp(measles_cases,
          units = "town",
          times = "time",
          t0 = min(measles_cases$time)-1/52,
          unit_statenames = unit_statenames,
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

measles_params <- rep(0,length=length(measles_paramnames))
names(measles_params) <- measles_paramnames

for(p in fixedParNames) measles_params[paste0(p,1)] <- basic_params[p]
for(p in expandedParNames) measles_params[paste0(p,1:U)] <- basic_params[p]

coef(m1) <- measles_params
m1
}
