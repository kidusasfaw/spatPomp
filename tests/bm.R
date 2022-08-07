library(spatPomp)
set.seed(2)

b_U <- 3
b_N <- 3

b_model <- bm(U=b_U,N=b_N) 

b_pf <- pfilter(b_model,Np=10)
logLik(b_pf)

extended <- TRUE
if(extended){

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

b_igirf_lookahead <- 1
b_igirf_ninter <- 2
b_igirf_np <- 10
b_igirf_nguide <- 10
b_igirf_ngirf <- 10

b_igirf_out <- igirf(b_model,
  Ngirf = b_igirf_ngirf,
  rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02,
    X1_0=ivp(0), X2_0=ivp(0), X3_0=ivp(0)
  ),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=b_igirf_np,
  Ninter = b_igirf_ninter,
  lookahead = b_igirf_lookahead,
  Nguide = b_igirf_nguide,
  kind = 'moment',
  verbose = FALSE
)

logLik(b_igirf_out)
}

