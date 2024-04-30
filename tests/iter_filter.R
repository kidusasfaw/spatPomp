
## mostly tested in various calls to iterated filter algorithms in bm.R

## here, we just test edge cases

library(spatPomp)
try(spatPomp:::perturbn.kernel.sd(rw.sd=rw_sd(rho=0.02,X1_0=ivp(0.02)),1:2,
  paramnames="JUNK"))

try(spatPomp:::perturbn.kernel.sd("JUNK"))

try(spatPomp:::perturbn.kernel.sd(rw.sd=rw_sd(0.02)))

try(spatPomp:::perturbn.kernel.sd(rw.sd=rw_sd(rho=1:10,X1_0=ivp(0.02)),1:2,
  paramnames=c("rho","X1_0")))

spatPomp:::safecall()

spatPomp:::perturbn.kernel.sd(
  rw.sd=matrix(c(0.01,0.02),nrow=2,ncol=2,
    dimnames=list(c("rho","X1_0"),NULL)), 1:2, paramnames=c("rho","X1_0"))

# test fraction > 1
frac_test <- spatPomp:::mif2.cooling("hyperbolic",fraction=1.5,ntimes=5)
frac_test(5,1)







