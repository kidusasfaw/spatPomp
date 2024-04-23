

library(spatPomp)
set.seed(1)

## a default measles model with fixed shared initial value parameters (IVPs)
## compiled into the object

m_U <- 2 
m_model <- measles(U=m_U)

range(time(m_model))

m <- window(m_model,1950,1955)

m_params <- c(R0=56.8,amplitude=0.554,gamma=30.4,sigma=28.9,mu=0.02,sigmaSE=0.02,rho=0.488,psi=0.116,g=100,alpha=1,iota=0,cohort=0.1)

m_pf <- pfilter(m_model,Np=20,params=m_params)
logLik(m_pf)

## setting IVPs

m2_model <- measles(U=m_U,fixed_ivps=FALSE)
m2_model <- window(m2_model,1950,1955)

m2_shared_params <- c(R0=56.8,amplitude=0.554,gamma=30.4,sigma=28.9,mu=0.02,sigmaSE=0.02,rho=0.488,psi=0.116,g=100,alpha=1,iota=0,cohort=0.1)
m2_unit_ivp_names <- c('S','E','I')
m2_statenames <- paste0(rep(m2_unit_ivp_names,each=m_U),1:m_U)
m2_ivp_names <- paste0(m2_statenames,"_0")
m2_ivps <- rep(c(0.032,0.00005,0.00004),each=m_U)
names(m2_ivps) <- m2_ivp_names
m2_params <- c(m2_shared_params,m2_ivps)

set.seed(1)
m2_pf <- pfilter(m2_model,Np=20,params=m2_params)
logLik(m2_pf)

##
## Note: the measles skeleton is correct only when there is no cohort effect,
## i.e., cohort=0.
##

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

m3_params <- m_params
m3_params["cohort"] <- 0
m3_igirf_lookahead <- 1
m3_igirf_ninter <- 2
m3_igirf_np <- 5
m3_igirf_nguide <- 5
m3_igirf_ngirf <- 2

m3_igirf_out <- igirf(m_model, Ngirf = m3_igirf_ngirf,
  params=m3_params,
  rw.sd=rw_sd(g=0.02),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=m3_igirf_np,
  Ninter = m3_igirf_ninter,
  lookahead = m3_igirf_lookahead,
  Nguide = m3_igirf_nguide,
  kind = 'moment',
  verbose = FALSE
)
logLik(m3_igirf_out)

## test error message
try(measles(U=1000))


## --------------------------------------------------------------------
## using measles to test spatPomp() covariate replacement functionality
## ____________________________________________________________________

m_covrep <- spatPomp(m,
  covar=spatPomp::measlesUK[
    measlesUK$city %in% m@unit_names,
    c("year","city","pop","births")]
)
for(slt in slotNames(m)) if(!identical(slot(m,slt),slot(m_covrep,slt))) print(slt)
## covariates are not the same since this replacement does not deal with
## birth delay