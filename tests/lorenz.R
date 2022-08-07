
library(spatPomp)
set.seed(3)

l_U <- 5
l_N <- 3

l <- lorenz(U=l_U,N=l_N)

l_pf <- pfilter(l,Np=10)
logLik(l_pf)

extended <- TRUE
if(extended){

## A call to igirf using the moment-based guide function can test compiled code for eunit_measure, munit_measure, vunit_measure, dunit_measure, runit_measure, rprocess, skeleton, rinit and partrans. 

l_igirf_lookahead <- 1
l_igirf_ninter <- 2
l_igirf_np <- 10
l_igirf_nguide <- 10
l_igirf_ngirf <- 10

l_igirf_out <- igirf(l,
  Ngirf = l_igirf_ngirf,
  rw.sd = rw.sd(F=0.02, sigma=0.02, tau=0.02,
    X1_0=ivp(0),X2_0=ivp(0), X3_0=ivp(0),X4_0=ivp(0)
  ),
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=l_igirf_np,
  Ninter = l_igirf_ninter,
  lookahead = l_igirf_lookahead,
  Nguide = l_igirf_nguide,
  kind = 'moment',
  verbose = FALSE
)

logLik(l_igirf_out)
}


