## extra tests of correctness requiring addition Monte Carlo intensity

set.seed(42)
library(spatPomp)


  print("Test he10 with towns_selected argument")
  h4 <- he10(U=4,towns_selected=c(1,2,11,12),
  basic_params = c(
    alpha =0.99,      iota=0,          R0=30,
    cohort=0.5,  amplitude=0.3,     gamma=52,
    sigma=52,           mu=0.02,  sigmaSE=0.05,
    rho=0.5,           psi=0.1,         g=800,
    S_0=0.036,         E_0=0.00007,   I_0=0.00006
  ))
  s4 <- simulate(h4,seed=27)
  obs(s4)[,1:2]
 