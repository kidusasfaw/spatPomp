library(spatPomp)
set.seed(2)

gbm_model <- gbm(U=3,N=2) 

gbm_pf <- pfilter(gbm_model,Np=10)
logLik(gbm_pf)

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

gbm_igirf_out <- igirf(gbm_model,
  Ngirf = 2,
  rw.sd = rw_sd(rho=0.02, sigma=0.02, tau=0.02),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=10,
  Ninter = 2,
  lookahead = 1,
  Nguide = 5,
  kind = 'moment',
  verbose = FALSE
)

logLik(gbm_igirf_out)


