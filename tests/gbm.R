library(spatPomp)
set.seed(2)

gbm_U <- 3
gbm_N <- 2

gbm_model <- bm(U=gbm_U,N=gbm_N) 

gbm_pf <- pfilter(gbm_model,Np=10)
logLik(gbm_pf)

extended <- TRUE
if(extended){

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

gbm_igirf_lookahead <- 1
gbm_igirf_ninter <- 2
gbm_igirf_np <- 10
gbm_igirf_nguide <- 5
gbm_igirf_ngirf <- 2

gbm_igirf_out <- igirf(gbm_model,
  Ngirf = gbm_igirf_ngirf,
  rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=gbm_igirf_np,
  Ninter = gbm_igirf_ninter,
  lookahead = gbm_igirf_lookahead,
  Nguide = gbm_igirf_nguide,
  kind = 'moment',
  verbose = FALSE
)

logLik(gbm_igirf_out)
}

