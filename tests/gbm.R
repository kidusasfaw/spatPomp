library(spatPomp)
set.seed(2)

gbm_model <- gbm(U=3,N=2) 

gbm_pf <- pfilter(gbm_model,Np=10)

paste("gbm pfilter loglik:",round(logLik(gbm_pf),10))

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

paste("gbm igirf loglik:", round(logLik(gbm_igirf_out),10))

## --------------------------------------------
## using gbm to test arma_benchmark()
## ____________________________________________


# unit test commented out for 0.33.0 since CRAN additional checks found
# that stats::arima() has reproducibility problems on an M1 mac.

# a1 <- arma_benchmark(gbm_model,order=c(1,0,0))

# paste("ARMA benchmark:", round(a1$total,10))




