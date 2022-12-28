
library(spatPomp)
set.seed(3)

l <- lorenz(U=5,N=3)

l_pf <- pfilter(l,Np=10)
paste("lorenz pfilter loglik: ", logLik(l_pf))

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

l_igirf <- igirf(l,
  Ngirf = 2,
  rw.sd = rw_sd(F=0.02, tau=0.02,X1_0=ivp(0),X2_0=ivp(0)),
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  Np=10,
  Ninter = 2,
  lookahead = 1,
  Nguide = 10,
  kind = 'bootstrap',
  verbose = FALSE
)
paste("lorenz igirf bootstrap hyperbolic loglik:", round(logLik(l_igirf),10))



