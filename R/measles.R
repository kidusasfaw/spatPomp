#' Measles in UK spatPomp generator
#'
#' Generate a spatPomp object for measles in the top-\code{U} most populous cities in England.
#' Model adapted from He et al. (2010) with gravity transport following Park and Ionides (2019).
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatPomp object.
#' @return A spatPomp object.
#' @examples
#' measles(7)
#' @export

measles <- function(U=6){

birth_lag <- 3*26  # delay until births hit susceptibles, in biweeks

# pre-vaccine biweekly measles reports for the largest 40 UK cities, sorted by size
data(measlesUK)
measlesUK$city<-as.character(measlesUK$city)

# Note: check for outliers, c.f. He et al (2010)


######## code for data cleaning: only re-run if dataset changes ######
if(0){
# datafile for measles spatPomp
# derived from measlesUKUS.csv from
# https://datadryad.org/resource/doi:10.5061/dryad.r4q34
# US data come from Project Tycho.
# England and Wales data are the city of London plus the largest 39 cities that were more than 50km from London.
# cases is reported measles cases per biweek
# births is estimated recruitment of susceptibles per biweek
  library(magrittr)
  library(dplyr)
  read.csv("../../measles/measlesUKUS.csv",stringsAsFactors=FALSE) %>% subset(country=="UK") -> x
  library(dplyr)
x %>%
  group_by(loc) %>%
  mutate(meanPop = mean(pop)) %>%
  ungroup() %>%
  arrange(desc(meanPop),decimalYear) -> x1
x1 %>% transmute(year=decimalYear,city=loc,cases=cases,pop=pop,births=rec) -> x2
  # the R package csv format
  # from https://cran.r-project.org/doc/manuals/R-exts.html#Data-in-packages
  write.table(file="measlesUK.csv",sep = ";",row.names=F,x2)
y <- x1[x1$decimalYear==1944,c("loc","lon","lat","meanPop")]
y1 <- transmute(y,city=loc,lon,lat,meanPop)
write.table(file="city_data_UK.csv",sep=";",row.names=F,y1)
}
####################################################################

cities <- unique(measlesUK$city)[1:U]
measles_cases <- measlesUK[measlesUK$city %in% cities,c("year","city","cases")]
measles_cases <- measles_cases[measles_cases$year>1949.99,]
measles_covar <- measlesUK[measlesUK$city %in% cities,c("year","city","pop","births")]
u <- split(measles_covar$births,measles_covar$city)
v <- sapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
measles_covar$lag_birthrate <- as.vector(v[,cities])*26
measles_covar$births<- NULL
measles_covarnames <- paste0(rep(c("pop","lag_birthrate"),each=U),1:U)

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

measles_globals <- Csnippet(
  paste0("const int U = ",U,"; \n ", v_by_g_C)
)

measles_unit_statenames <- c('S','E','I','R','C','W')
#measles_unit_statenames <- c('S','E','I','R','Acc','C','W')

measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
measles_IVPnames <- paste0(measles_statenames[1:(4*U)],"_0")
measles_RPnames <- c("alpha","iota","R0","cohort","amplitude","gamma","sigma","mu","sigmaSE","rho","psi","g")
measles_paramnames <- c(measles_RPnames,measles_IVPnames)

measles_rprocess <- Csnippet('
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
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

  // pre-computing this saves substantial time
  for (u = 0 ; u < U ; u++) {
    powVec[u] = pow(I[u]/pop[u],alpha);
  }

  for (u = 0 ; u < U ; u++) {

    // cohort effect
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
      br = cohort*lag_birthrate[u]/dt + (1-cohort)*lag_birthrate[u];
    else
      br = (1.0-cohort)*lag_birthrate[u];

    // expected force of infection
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

    // These rates could be outside the u loop if all parameters are shared between units
    rate[1] = mu;			    // natural S death
    rate[2] = sigma;		  // rate of ending of latent stage
    rate[3] = mu;			    // natural E death
    rate[4] = gamma;		  // recovery
    rate[5] = mu;			    // natural I death

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
    W[u] += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
    C[u] += trans[4];           // true incidence

    // For girf, we want to override the default behavior of the rprocess
    // function, which sets accumulator variables (in this case C) to 0
    // at each call. We only want this if the call is during an observation
    // time.
    //if(obstime){
       //Acc[u] = trans[4];
       //C[u] = Acc[u];
    //} else{
       //Acc[u] += trans[4];
       //C[u] = Acc[u];
    //}
  }
')


measles_dmeasure <- Csnippet("
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

measles_rmeasure <- Csnippet("
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

measles_unit_dmeasure <- Csnippet('
                       double m = rho*C;
                       double v = m*(1.0-rho+psi*psi*m);
                       double tol = 1.0e-18;
                       if (cases > 0.0) {
                         lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                       } else {
                           lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                       }
                       ')

measles_unit_emeasure <- Csnippet("
ey = rho*C;
")

measles_unit_vmeasure <- Csnippet("
double m;
m = rho*C;
vc = m*(1.0-rho+psi*psi*m);
")

measles_unit_mmeasure <- Csnippet("
double binomial_var;
double m;
m = rho*C;
binomial_var = rho*(1-rho)*C;
if(vc > binomial_var) M_psi = sqrt(vc - binomial_var)/m;
")

measles_rinit <- Csnippet("
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
  //double *Acc = &Acc1;
  const double *S_0 = &S1_0;
  const double *E_0 = &E1_0;
  const double *I_0 = &I1_0;
  const double *R_0 = &R1_0;
  //const double *Acc_0 = &Acc1_0;
  const double *pop = &pop1;
  double m;
  int u;
  for (u = 0; u < U; u++) {
    m = pop[u]/(S_0[u]+E_0[u]+I_0[u]+R_0[u]);
    S[u] = nearbyint(m*S_0[u]);
    E[u] = nearbyint(m*E_0[u]);
    I[u] = nearbyint(m*I_0[u]);
    R[u] = nearbyint(m*R_0[u]);
    W[u] = 0;
    C[u] = 0;
    //Acc[u] = Acc_0[u];
  }
")

measles_skel <- Csnippet('
  double beta, br, seas, foi;
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
  double *DS = &DS1;
  double *DE = &DE1;
  double *DI = &DI1;
  double *DR = &DR1;
  //double *DAcc = &DAcc1;
  double *DC = &DC1;
  double *DW = &DW1;
  double powVec[U];
  //double *Acc = &Acc1;
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
    //if(obstime != 1){
       //C[u] = Acc[u];
    //}
    // cannot readily put the cohort effect into a vectorfield for the skeleton
    // therefore, we ignore it here.
    // this is okay as long as the skeleton is being used for short-term forecasts
    br = lag_birthrate[u];

    foi = pow( (I[u]+iota)/pop[u],alpha);
    for (v=0; v < U ; v++) {
      if(v != u)
        foi += g * v_by_g[u][v] * (powVec[v] - powVec[u]) / pop[u];
    }

    DS[u] = br - (beta*foi + mu)*S[u];
    DE[u] = beta*foi*S[u] - (sigma+mu)*E[u];
    DI[u] = sigma*E[u] - (gamma+mu)*I[u];
    DR[u] = gamma*I[u] - mu*R[u];
    DW[u] = 0;
    DC[u] = gamma*I[u];
    //DAcc[u] = 0;
  }
')


spatPomp(measles_cases,
        units = "city",
        times = "year",
        t0 = min(measles_cases$year)-1/26,
        unit_statenames = measles_unit_statenames,
        covar = measles_covar,
        tcovar = "year",
        rprocess=euler(measles_rprocess, delta.t=2/365),
        skeleton=vectorfield(measles_skel),
        accumvars = c(paste0("C",1:U),paste0("W",1:U)),
        paramnames=measles_paramnames,
        covarnames=measles_covarnames,
        globals=measles_globals,
        rinit=measles_rinit,
        dmeasure=measles_dmeasure,
        unit_emeasure=measles_unit_emeasure,
        unit_mmeasure=measles_unit_mmeasure,
        unit_vmeasure=measles_unit_vmeasure,
        rmeasure=measles_rmeasure,
        unit_dmeasure=measles_unit_dmeasure
)
}

#' @export
girfd_measles <- function(U=5, N = 10, Np = 100, Nguide = 50){
  # Get real data for U=40 so we can simulate using the resulting spatPomp
  measles_sim_U <- 40
  measles_uk <- measles(measles_sim_U)
  read.csv(text="
,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
LONDON,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
BIRMINGHAM,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
LIVERPOOL,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
MANCHESTER,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
LEEDS,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
SHEFFIELD,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
BRISTOL,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
NOTTINGHAM,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
HULL,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
BRADFORD,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
",stringsAsFactors=FALSE,row.names=1) -> he10_mles

  # Set the parameters for the simulation
  measles_unit_statenames <- c('S','E','I','R', 'C','W')
  measles_statenames <- paste0(rep(measles_unit_statenames,each=measles_sim_U),1:measles_sim_U)
  measles_IVPnames <- paste0(measles_statenames[1:(4*measles_sim_U)],"_0")
  measles_RPnames <- c("alpha","iota","R0","cohort","amplitude",
                       "gamma","sigma","mu","sigmaSE","rho","psi","g")
  measles_paramnames <- c(measles_RPnames,measles_IVPnames)
  measles_params <- rep(NA,length(measles_paramnames))
  names(measles_params) <- measles_paramnames
  city_params <- unlist(he10_mles["LONDON",])
  measles_params[measles_RPnames] <- c(city_params,g=100)[measles_RPnames]
  measles_params[paste0("S",1:measles_sim_U,"_0")] <-city_params["S_0"]
  measles_params[paste0("E",1:measles_sim_U,"_0")] <-city_params["E_0"]
  measles_params[paste0("I",1:measles_sim_U,"_0")] <-city_params["I_0"]
  measles_params[paste0("R",1:measles_sim_U,"_0")] <-city_params["R_0"]
  # measles_params[paste0("Acc",1:measles_sim_U,"_0")] <- 0

  # Perform a 40-city simulation which will then be subsetted
  measles_sim <- simulate(measles_uk,params=measles_params)

  # Subsetting function
  measles_subset <- function(m_U,m_N){
    m <- measles(U=m_U)
    m@data <- measles_sim@data[1:m_U,1:m_N]
    time(m) <- measles_sim@times[1:m_N]
    m_statenames <- paste0(rep(measles_unit_statenames,each=m_U),1:m_U)
    m_IVPnames <- paste0(m_statenames[1:(4*m_U)],"_0")
    m_paramnames <- c(measles_RPnames,m_IVPnames)
    m_params <- measles_params[names(measles_params)%in%m_paramnames]
    coef(m) <- m_params
    return(m)
  }

  # Get simulated data for U cities and N times
  m <- measles_subset(m_U=U, m_N=N)

  # gird_spatPomp object creation requirements
  measles_Ninter <- length(spat_units(m))
  measles_lookahead <- 1
  measles_Nguide <- Nguide
  measles_Np <- Np
  measles_tol <- 1e-300
  # measles_h <- function(state.vec, param.vec){
  #   ix<-grep("C",names(state.vec))
  #   param.vec["rho"]*state.vec[ix]
  # }
  # measles_theta.to.v <- function(meas.mean, param.vec){
  #   meas.mean*(1-param.vec["rho"]) + (meas.mean*param.vec["psi"])^2
  # }
  # measles_v.to.theta <- function(var, state.vec, param.vec){
  #   cases <- state.vec[grep("C",names(state.vec))]
  #   binomial_var <- param.vec['rho']*(1-param.vec["rho"])*cases
  #   param_fit <- param.vec
  #   if(var > binomial_var) {
  #     param_fit["psi"] <- sqrt(var - binomial_var) / (param.vec["rho"]*cases)
  #   }
  #   param_fit
  # }

  # Output girfd_spatPomp object
  new(
    "girfd_spatPomp",
    m,
    Ninter=measles_Ninter,
    Nguide=measles_Nguide,
    lookahead=measles_lookahead,
    cond.loglik = array(data=numeric(0),dim=c(0,0)),
    Np = as.integer(measles_Np),
    tol= measles_tol,
    loglik=as.double(NA)
  )

}

#' @export
asifird_measles <- function(U=5,
                            N = 10,
                            islands = 50,
                            Np = 10,
                            nbhd = function(object, time, unit) {
                              nbhd_list <- list()
                              if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
                              if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
                              return(nbhd_list)
                            },
                            Ninter = U
                            ){
  # Get real data for U=40 so we can simulate using the resulting spatPomp
  measles_sim_U <- 40
  measles_uk <- measles(measles_sim_U)
  read.csv(text="
,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
LONDON,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
BIRMINGHAM,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
LIVERPOOL,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
MANCHESTER,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
LEEDS,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
SHEFFIELD,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
BRISTOL,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
NOTTINGHAM,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
HULL,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
BRADFORD,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
",stringsAsFactors=FALSE,row.names=1) -> he10_mles

  # Set the parameters for the simulation
  measles_unit_statenames <- c('S','E','I','R', 'C','W')
  measles_statenames <- paste0(rep(measles_unit_statenames,each=measles_sim_U),1:measles_sim_U)
  measles_IVPnames <- paste0(measles_statenames[1:(4*measles_sim_U)],"_0")
  measles_RPnames <- c("alpha","iota","R0","cohort","amplitude",
                       "gamma","sigma","mu","sigmaSE","rho","psi","g")
  measles_paramnames <- c(measles_RPnames,measles_IVPnames)
  measles_params <- rep(NA,length(measles_paramnames))
  names(measles_params) <- measles_paramnames
  city_params <- unlist(he10_mles["LONDON",])
  measles_params[measles_RPnames] <- c(city_params,g=100)[measles_RPnames]
  measles_params[paste0("S",1:measles_sim_U,"_0")] <-city_params["S_0"]
  measles_params[paste0("E",1:measles_sim_U,"_0")] <-city_params["E_0"]
  measles_params[paste0("I",1:measles_sim_U,"_0")] <-city_params["I_0"]
  measles_params[paste0("R",1:measles_sim_U,"_0")] <-city_params["R_0"]
  # measles_params[paste0("Acc",1:measles_sim_U,"_0")] <- 0

  # Perform a 40-city simulation which will then be subsetted
  measles_sim <- simulate(measles_uk,params=measles_params)

  # Subsetting function
  measles_subset <- function(m_U,m_N){
    m <- measles(U=m_U)
    m@data <- measles_sim@data[1:m_U,1:m_N]
    time(m) <- measles_sim@times[1:m_N]
    m_statenames <- paste0(rep(measles_unit_statenames,each=m_U),1:m_U)
    m_IVPnames <- paste0(m_statenames[1:(4*m_U)],"_0")
    m_paramnames <- c(measles_RPnames,m_IVPnames)
    m_params <- measles_params[names(measles_params)%in%m_paramnames]
    coef(m) <- m_params
    return(m)
  }

  # Get simulated data for U cities and N times
  m <- measles_subset(m_U=U, m_N=N)

  # gird_spatPomp object creation requirements
  measles_Ninter <- Ninter
  measles_Np <- Np
  measles_tol <- 1e-300
  measles_nbhd <- nbhd
  measles_islands <- islands

  # Output girfd_spatPomp object
  new(
    "asifird_spatPomp",
    m,
    Np = as.integer(measles_Np),
    Ninter=as.integer(measles_Ninter),
    islands=as.integer(measles_islands),
    nbhd=measles_nbhd,
    tol= measles_tol,
    loglik=as.double(NA)
  )

}

