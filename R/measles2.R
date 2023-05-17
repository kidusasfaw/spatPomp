#' Measles in UK: spatPomp generator with shared or unit-specific parameters
#'
#' Generate a spatPomp object for measles in the top-\code{U} most populous cities in England and Wales.
#' The model is adapted from He et al. (2010) with gravity transport following Park and Ionides (2019).
#' The structure of this spatPomp is designed to accommodate shared and unit-specific parameters.
#' If carrying out spatiotemporal iterated filtering for shared parameters via ibpf, it is necessary to
#' have a unit-specific expansion and so these parameters should be included in expandedParNames.
#' This model and data correspond to the biweekly analysis of Park and Ionides (2020) and
#' Ionides et al (2021). There are small differences with the weekly model and data of
#' He et al (2010) and Ionides, Ning and Wheeler (2022).
#'
#' @param U An integer from 1 to 40 specifying the number of cities to be represented in the spatPomp object.
#' @importFrom utils data read.csv write.table
#' @param N An integer from 1 to 391 specifying the number of time points.
#' @param basic_params A named vector used to specify shared parameters or unit-specific parameters
#' having common values for each unit.
#' @param dt a numeric (in unit of years) that is used as the Euler time-increment for simulating measles data
#' @param expandedParNames specifies parameters that are defined for each unit. This also allows unit perturbations for a parameter with a value shared across units.
#' @param contractedParNames specifies parameters having a shared value across units. Remaining parameters that are neither expanded nor contracted are considered fixed, and will not have a transformation defined for them.
#' @param simulated determines whether to return a simulation from the model or the
#' UK measles data
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled measles POMP model.
#' @references
#'
#' \he2010
#'
#' \park2020
#'
#' \ionides2021
#'
#' \ionides2022
#'
#' @examples
#' # Complete examples are provided in the package tests
#' \dontrun{
#' m <- measles2(U = 5)
#' # See all the model specifications of the object
#' spy(m)
#' }
#' @export

measles2 <- function(U=6,dt=2/365,N=391,
  expandedParNames=c("R0", "c", "A", "muIR",
    "muEI", "sigmaSE", "rho", "psi", "g", "S_0", "E_0", "I_0"),
  contractedParNames=NULL,
  simulated=FALSE,
  basic_params =c(
    alpha = 0.98, # exponent on infectives
    iota = 0.1, # immigration (number of infectives)
    R0 = 30,
    c = 0.3,    # cohort fraction
    A = 0.5,    # amplitude of seasonal variation in transmission rate
    muIR = 52,  # recovery rate (year^-1)
    muEI = 52,  # inverse of mean latent duration (year)
    muD = 0.02, # inverse of life expectancy (year)
    sigmaSE = 0.15, # process noise (year^(1/2))
    rho = 0.5,  # reporting rate
    psi = 0.15, # measurement overdispersion
    g = 400,    # movement model gravitational constant
    S_0 = 0.032,   # initial fraction susceptible
    E_0 = 0.00005, 
    I_0 = 0.00004
  )
){
  if(U>40) stop("U <= 40")
  if(N>391) stop("N <= 391")
  birth_lag <- 4*26  # delay until births hit susceptibles, in biweeks

  # pre-vaccine biweekly measles reports for the largest 40 UK cities, sorted by size
  measlesUK <- spatPomp::measlesUK
  measlesUK$city<-as.character(measlesUK$city)
  city_data_UK <- spatPomp::city_data_UK

  cities <- unique(measlesUK$city)[1:U]
  selected <- (measlesUK$city %in% cities) & (measlesUK$year>1949.99) &
    (measlesUK$year<1950.01+(N-1)/26)    
  measles_cases <- measlesUK[selected,c("year","city","cases")]
  covar_selected <- (measlesUK$city %in% cities) & (measlesUK$year<1950.01+(N-1)/26) 
  measles_covar <- measlesUK[covar_selected,c("year","city","pop","births")]
  u <- split(measles_covar$births,measles_covar$city)
  v <- sapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
  measles_covar$lag_birthrate <- as.vector(v[,cities])*26
  measles_covar$births <- NULL
  names(measles_covar) <- c("year","city","P","lag_birthrate")

  # Haversine formula for great circle distance between two points
  # on a sphere radius r. Here, r defaults to a mean radius for the earth, in miles.
  distGreatCircle <- function(p1, p2, r = 3963.191) {
    Lon1 <- p1[,1]*pi/180
    Lat1 <- p1[,2]*pi/180
    Lon2 <- p2[,1]*pi/180
    Lat2 <- p2[,2]*pi/180
    a <- sin((Lat2-Lat1)/2)^2 + cos(Lat1)*cos(Lat2)*sin((Lon2-Lon1)/2)^2
    atan2(sqrt(a), sqrt(1 - a)) * 2 * r
  }

  lon_lat <- city_data_UK[1:U,c("lon","lat")]
  dmat <- matrix(0,U,U)
  for(u1 in 1:U) {
    for(u2 in 1:U) {
      dmat[u1,u2] <- round(distGreatCircle(lon_lat[u1,],lon_lat[u2,]),1)
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

  basic_RPnames <- c("alpha","iota","psi","R0","muIR","muEI","sigmaSE","c","A","muD","rho","g")
  basic_IVPnames <- c("S_0", "E_0", "I_0")
  basicParNames <- c(basic_RPnames,basic_IVPnames)

  fixedParNames <- setdiff(basicParNames,c(expandedParNames,contractedParNames))

  set_unit_specific <- Csnippet(paste0("const int ", expandedParNames,
    "_unit = 1;\n", collapse=" "))
  set_shared <- Csnippet(paste0("const int ", c(contractedParNames,fixedParNames),
    "_unit = 0;\n", collapse=" "))

  measles_globals <- Csnippet(
    paste(v_by_g_C, set_unit_specific, set_shared, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  measles_paramnames <- unlist(c(
    lapply(expandedParNames, function(x,U) paste0(x,1:U),U),
    lapply(c(contractedParNames,fixedParNames),function(x) paste0(x,1))
  ))

  measles_rprocess <- spatPomp_Csnippet(
    unit_statenames = c('S','E','I','C'),
    unit_covarnames = c('P','lag_birthrate'),
    unit_paramnames = c('alpha','iota','R0','c','A','muIR','muEI','muD','sigmaSE','g'),
    code="
      int BS=0, SE=1, SD=2, EI=3, ED=4, IR=5, ID=6;
      double seas, Ifrac, dw;
      double mu[7], dN[7];
      double powVec[U];
      int u,v;
      double day = (t-floor(t))*365;

      for (u = 0 ; u < U ; u++) {

        // needed for the Ensemble Kalman filter
        // or other methods making real-valued perturbations to the state
        // reulermultinom requires integer-valued double type for states
        S[u] = S[u]>0 ? floor(S[u]) : 0;
        E[u] = E[u]>0 ? floor(E[u]) : 0;
        I[u] = I[u]>0 ? floor(I[u]) : 0;

        // pre-computing this saves substantial time
        powVec[u] = pow(I[u]/P[u],alpha[u*alpha_unit]);
      }

      for (u = 0 ; u < U ; u++) {

        seas = (day>=7&&day<=100)||(day>=115&&day<=199)||(day>=252&&day<=300)||(day>=308&&day<=356)
          ? 1.0 + A[u*A_unit] * 0.2411/0.7589 : 1.0 - A[u*A_unit];

        // cohort effect
        if (fabs(day-251.0) < 0.5*dt*365)
          mu[BS] = c[u*c_unit]*lag_birthrate[u]/dt + (1-c[u*c_unit])*lag_birthrate[u];
        else
          mu[BS] = (1.0-c[u*c_unit])*lag_birthrate[u];

        // we follow Park and Ionides (2019) and raise pop to the alpha power
        // He et al (2010) did not do this.
        Ifrac = pow( (I[u]+iota[u*iota_unit])/P[u],alpha[u*alpha_unit]);
        for (v=0; v < U ; v++) {
           if(v != u)
            Ifrac += g[u*g_unit] * v_by_g[u][v] * (powVec[v] - powVec[u]) / P[u];
        }

        dw = rgammawn(sigmaSE[u*sigmaSE_unit],dt);
        mu[SE] = R0[u*R0_unit]*(muIR[u*muIR_unit]+muD[u*muD_unit])*seas*Ifrac*dw/dt;
        mu[SD] = muD[u*muD_unit];		
        mu[EI] = muEI[u*muEI_unit];	
        mu[ED] = muD[u*muD_unit];		
        mu[IR] = muIR[u*muIR_unit];	
        mu[ID] = muD[u*muD_unit];		
 
        // transitions between classes
        reulermultinom(2,S[u],&mu[SE],dt,&dN[SE]);  // SE and SD transitions
        reulermultinom(2,E[u],&mu[EI],dt,&dN[EI]);  // EI and ED transitions
        reulermultinom(2,I[u],&mu[IR],dt,&dN[IR]);  // IR and ID transitions
        dN[BS] = rpois(mu[BS]*dt);
 
        S[u] += dN[BS]  - dN[SE] - dN[SD];
        E[u] += dN[SE] - dN[EI] - dN[ED];
        I[u] += dN[EI] - dN[IR] - dN[ID];
        C[u] += dN[IR];          
      }
    "
  )

  measles_dmeasure <- spatPomp_Csnippet(
    unit_statenames = 'C',
    unit_obsnames = 'cases',
    unit_paramnames = c('rho','psi'),
    code="
      double m,v;
      double tol = 1e-300;
      double mytol = 1e-5;
      int u;

      lik = 0;
      for (u = 0; u < U; u++) {
        m = rho[u*rho_unit]*(C[u]+mytol);
        v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
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
    "
  )

  measles_rmeasure <- spatPomp_Csnippet(
    method='rmeasure',
    unit_paramnames=c('rho','psi'),
    unit_statenames='C',
    unit_obsnames='cases',
    code="
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
    "
  )

  measles_dunit_measure <- spatPomp_Csnippet(
    unit_paramnames=c('rho','psi'),
    code="
      double mytol = 1e-5;
      double m = rho[u*rho_unit]*(C+mytol);
      double v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
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
    "
  )

  measles_eunit_measure <- spatPomp_Csnippet(
    unit_paramnames='rho',
    code="
      ey = rho[u*rho_unit]*C;
    "
  )

  measles_vunit_measure <- spatPomp_Csnippet(
    unit_paramnames=c('rho','psi'),
    code="
      double mytol = 1e-5;
      double m;
      m = rho[u*rho_unit]*(C+mytol);
      vc = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
    "
  )

  measles_rinit <- spatPomp_Csnippet(
    unit_paramnames = c('S_0','E_0','I_0'),
    unit_statenames = c('S','E','I','C'),
    unit_covarnames = 'P',
    code="
      double m;
      int u;
      for(u = 0; u < U; u++) {
        m = (float)(P[u]);
        S[u] = nearbyint(m*S_0[u*S_0_unit]);
        I[u] = nearbyint(m*I_0[u*I_0_unit]);
        E[u] = nearbyint(m*E_0[u*E_0_unit]);
        C[u] = 0;
      }
    "
  )

  logParNames <- c("muEI", "muIR", "sigmaSE", "psi", "R0", "g","iota","muD")
  logExpandedParNames <- intersect(logParNames,expandedParNames)
  logContractedParNames <- intersect(logParNames,contractedParNames)

  logitParNames <- c("A", "alpha","c","rho","S_0", "E_0", "I_0")
  logitExpandedParNames <- intersect(logitParNames,expandedParNames)
  logitContractedParNames <- intersect(logitParNames,contractedParNames)

  log_names <- c(
    unlist(lapply(logExpandedParNames, function(x,U) paste0(x,1:U),U)),
    unlist(lapply(logContractedParNames, function(x) paste0(x,1)))
  )

  logit_names <- c(
    unlist(lapply(logitExpandedParNames, function(x,U) paste0(x,1:U),U)),
    unlist(lapply(logitContractedParNames, function(x) paste0(x,1)))
  )
  
  measles_partrans <- parameter_trans(log=log_names,logit=logit_names)

  m1 <-  spatPomp(measles_cases,
    units = "city",
    times = "year",
    t0 = min(measles_cases$year)-1/26,
    unit_statenames = c('S','E','I','C'),
    covar = measles_covar,
    rprocess=euler(measles_rprocess, delta.t=dt),
    unit_accumvars = 'C',
    paramnames=measles_paramnames,
    partrans=measles_partrans,
    globals=measles_globals,
    rinit=measles_rinit,
    dmeasure=measles_dmeasure,
    rmeasure=measles_rmeasure,
    dunit_measure=measles_dunit_measure,
    eunit_measure=measles_eunit_measure,    
    vunit_measure=measles_vunit_measure           
  )

  measles_params <- rep(0,length=length(measles_paramnames))
  names(measles_params) <- measles_paramnames

  for(p in c(fixedParNames,contractedParNames)) measles_params[paste0(p,1)] <- basic_params[p]
  for(p in expandedParNames) measles_params[paste0(p,1:U)] <- basic_params[p]

  coef(m1) <- measles_params
  if(simulated) simulate(m1) else m1
}
