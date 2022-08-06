

library(spatPomp)
set.seed(1)

## a default measles model with fixed shared initial value parameters (IVPs)
## compiled into the object

m_U <- 2 
m <- measles(U=m_U)

range(time(m))

m <- window(m,1950,1955)

m_params <- c(R0=56.8,amplitude=0.554,gamma=30.4,sigma=28.9,mu=0.02,sigmaSE=0.02,rho=0.488,psi=0.116,g=100,alpha=1,iota=0,cohort=0.1)

m_pf <- pfilter(m,Np=20,params=m_params)
logLik(m_pf)

## setting IVPs

m2 <- measles(U=m_U,fixed_ivps=FALSE)
m2 <- window(m2,1950,1955)

m2_shared_params <- c(R0=56.8,amplitude=0.554,gamma=30.4,sigma=28.9,mu=0.02,sigmaSE=0.02,rho=0.488,psi=0.116,g=100,alpha=1,iota=0,cohort=0.1)
m2_unit_ivp_names <- c('S','E','I')
m2_statenames <- paste0(rep(m2_unit_ivp_names,each=m_U),1:m_U)
m2_ivp_names <- paste0(m2_statenames,"_0")
m2_ivps <- rep(c(0.032,0.00005,0.00004),each=m_U)
names(m2_ivps) <- m2_ivp_names
m2_params <- c(m2_shared_params,m2_ivps)

set.seed(1)
m2_pf <- pfilter(m2,Np=20,params=m2_params)
logLik(m2_pf)

##
## Note: the measles skeleton is correct only when there is no cohort effect,
## i.e., cohort=0.
##



