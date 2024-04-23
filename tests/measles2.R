
library(spatPomp)

## test error messages
try(measles2(U=1000,N=5))
try(measles2(U=2,N=1000))

i <- 1
DEBUG=FALSE
U <- switch(i,4,10,40)
N <- switch(i,2,50,391)
Np <- switch(i,10,100,1000)

m1 <- measles2(U=U,N=N)
if(DEBUG){
  plot(simulate(m1),log=T)
}

# test for all parameters expanded, by default
set.seed(1)
s1 <- simulate(m1)
head(obs(s1))

if(DEBUG){
par1 <- coef(m1)
par1[paste0('iota',1:U)] <- 10
 plot(simulate(m1,params=par1),ty="l",log=T)
}

# test for all parameters contracted, e.g., for mif2 with
# all shared parameters
m2 <- measles2(U=U,N=N,expandedParNames=NULL,
  contractedParNames=c("R0", "c", "A", "muIR",
    "muEI", "sigmaSE", "rho", "psi", "g", "S_0", "E_0", "I_0")
)

set.seed(1)
s2 <- simulate(m2)

# test for all parameters contracted, e.g., for mif2 with
# all shared parameters
m2 <- measles2(U=U,N=N,expandedParNames=NULL,
  contractedParNames=c("R0", "c", "A", "muIR",
    "muEI", "sigmaSE", "rho", "psi", "g", "S_0", "E_0", "I_0")
)

# test for both expanded and contracted parameters
m3 <- measles2(U=U,N=N,expandedParNames=c("A","muIR"),
  contractedParNames=c("R0", "c", "S_0")
)

set.seed(1)
s3 <- simulate(m3)

if(any(obs(s2)!=obs(s1)))stop("s1 and s2 should be identical")
if(any(obs(s3)!=obs(s1)))stop("s1 and s3 should be identical")

partrans(m1,coef(m1),dir="toEst")
partrans(m2,coef(m2),dir="toEst")
partrans(m3,coef(m3),dir="toEst")

set.seed(2)
pf1 <- pfilter(m1,Np=Np)
logLik(pf1)
set.seed(2)
pf2 <- pfilter(m2,Np=Np)
logLik(pf2)
set.seed(2)
pf3 <- pfilter(m3,Np=Np)
logLik(pf3)

e1 <- enkf(m1,Np=Np)
logLik(e1)

if(DEBUG){
  # compare to measles()
  n1 <- measles(U=40)
  coef(n1) <- c(
    alpha = 1,
    iota = 0,  
    R0 = 30,
    cohort = 0,
    amplitude = 0.5,
    gamma = 52,
    sigma = 52,
    mu = 0.02,
    sigmaSE = 0.15, 
    rho = 0.5,
    psi = 0.15,
    g = 400,
    S_0 = 0.032, 
    E_0 = 0.00005, 
    I_0 = 0.00004
  )
  t1 <- simulate(n1)
  plot(simulate(n1),ty="l",log=T)
}

