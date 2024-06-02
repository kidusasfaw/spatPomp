library(spatPomp)
set.seed(2)

# Note: we need N >= 5 here to make the arma_benchmark stable
gbm_model <- gbm(U=2,N=5) 

gbm_pf <- pfilter(gbm_model,Np=5)

paste("gbm pfilter loglik:",round(logLik(gbm_pf),3))

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

gbm_igirf_out <- igirf(gbm_model,
  Ngirf = 1,
  rw.sd = rw_sd(rho=0.02, sigma=0.02, tau=0.02),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=3,
  Ninter = 1,
  lookahead = 1,
  Nguide = 4,
  kind = 'moment',
  verbose = FALSE
)

paste("gbm igirf loglik:", round(logLik(gbm_igirf_out),3))

## --------------------------------------------
## using gbm to test arma_benchmark()
## ____________________________________________

a1 <- arma_benchmark(gbm_model,order=c(1,0,0))

paste("ARMA benchmark:", round(a1$total,3))




