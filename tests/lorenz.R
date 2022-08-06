
library(spatPomp)
set.seed(2)

l_U <- 5
l_N <- 5

l <- lorenz(U=l_U,N=l_N)

l_pf <- pfilter(l,Np=10)
logLik(l_pf)



