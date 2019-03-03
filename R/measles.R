#' Measles in UK spatpomp generator
#'
#' Generate a spatpomp object for measles in the top-\code{U} most populous cities in England.
#' Model adapted from He et al. (2010) with gravity transport following Park and Ionides (2019).
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatpomp object.
#' @return A spatpomp object.
#' @examples
#' measles(7)
measles <- function(U=6){

birth_lag <- 3*26  # delay until births hit susceptibles, in biweeks

# pre-vaccine biweekly measles reports for the largest 40 UK cities, sorted by size
data(measlesUK)
measlesUK$city<-as.character(measlesUK$city)

# Note: check for outliers, c.f. He et al (2010)


######## code for data cleaning: only re-run if dataset changes ######
if(0){
# datafile for measles spatpomp
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
measles_statenames <- paste0(rep(measles_unit_statenames,each=U),1:U)
measles_IVPnames <- paste0(measles_statenames[1:(4*U)],"_0")
measles_RPnames <- c("alpha","iota","R0","cohort","amplitude","gamma","sigma","mu","sigmaSE","rho","psi","g")
measles_paramnames <- c(measles_RPnames,measles_IVPnames)

measles_rprocess <- Csnippet("
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
  const double *pop = &pop1;
  const double *lag_birthrate = &lag_birthrate1;
  int u,v;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;

  for (u = 0 ; u < U ; u++) {

    // cohort effect
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
      br = cohort*lag_birthrate[u]/dt + (1-cohort)*lag_birthrate[u];
    else
      br = (1.0-cohort)*lag_birthrate[u];

    // expected force of infection
    foi = pow( (I[u]+iota)/pop[u],alpha);
    // Do we still need iota in a spatPomp version?
    // See also discrepancy between Joonha and Daihai versions
    // Daihai didn't raise pop to the alpha power

    for (v=0; v < U ; v++) {
      if(v != u)
        foi += g * v_by_g[u][v] * (pow(I[v]/pop[v],alpha) - pow(I[u]/pop[u],alpha)) / pop[u];
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
  }
")

measles_rinit <- Csnippet("
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
  const double *S_0 = &S1_0;
  const double *E_0 = &E1_0;
  const double *I_0 = &I1_0;
  const double *R_0 = &R1_0;
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
  }
")

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

measles_rinit <- Csnippet("
  double *S = &S1;
  double *E = &E1;
  double *I = &I1;
  double *R = &R1;
  double *C = &C1;
  double *W = &W1;
  const double *S_0 = &S1_0;
  const double *E_0 = &E1_0;
  const double *I_0 = &I1_0;
  const double *R_0 = &R1_0;
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
  }
")

spatpomp(measles_cases,
                    units = "city",
                    times = "year",
                    t0 = min(measles_cases$year)-1/26,
                    unit_statenames = measles_unit_statenames,
                    covar = measles_covar,
                    tcovar = "year",
                    rprocess=euler(measles_rprocess, delta.t=2/365),
                    accumvars = c(paste0("C",1:U),paste0("W",1:U)),
                    paramnames=measles_paramnames,
                    covarnames=measles_covarnames,
                    globals=measles_globals,
                    rinit=measles_rinit,
                    dmeasure=measles_dmeasure,
                    rmeasure=measles_rmeasure,
                    unit_dmeasure=measles_unit_dmeasure
)

}
