library(spatPomp)
context("test methods on simple Brownian motion")

# create the BM object
U = 3; N = 10
bm3 <- bm(U = U, N = N)

# compute distance matrix to compute true log-likelihood
dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
dmat <- matrix(0,U,U)
for(u in 1:U) {
  for(v in 1:U) {
    dmat[u,v] <- dist(u,v)
  }
}

# compute the true log-likelihood
rootQ = coef(bm3)["rho"]^dmat * coef(bm3)["sigma"]
loglik.true <- pomp:::kalmanFilter(
  t=1:N,
  y=obs(bm3),
  X0=rinit(bm3),
  A= diag(length(spat_units(bm3))),
  Q=rootQ%*%rootQ,
  C=diag(1,nrow=nrow(dmat)),
  R=diag(coef(bm3)["tau"]^2, nrow=nrow(dmat))
)$loglik

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# h just extracts measured components
girf.h <- function(state.vec, param.vec){
  # find index matching unit_statename
  ix<-grep('X',names(state.vec))
  # no reporting ratio in model (equiv. to 1)
  state.vec[ix]
}

girf.theta.to.v <- function(meas.mean, param.vec){
  param.vec['tau']^2
}

girf.v.to.theta <- function(var, state.vec, param.vec){
  param.vec['tau'] <- sqrt(var)
  param.vec
}

girf.loglik <- replicate(10,logLik(girf(bm3,
                    Np = 100,
                    Ninter = length(spat_units(bm3)),
                    lookahead = 1,
                    Nguide = 50,
                    h = girf.h,
                    theta.to.v = girf.theta.to.v,
                    v.to.theta = girf.v.to.theta)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asif.nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

asif.loglik <- replicate(10,logLik(asif(bm3,
                           islands = 50,
                           Np = 20,
                           nbhd = asif.nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asifir.loglik <- replicate(10,logLik(asifir(bm3,
                        islands = 50,
                        Np=20,
                        nbhd = asif.nbhd,
                        Ninter = length(spat_units(bm3)),
                        h = girf.h,
                        theta.to.v = girf.theta.to.v,
                        v.to.theta = girf.v.to.theta)))

test_that("ASIF, ASIFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf.loglik) - loglik.true), 2)
  expect_lt(abs(logmeanexp(asif.loglik) - loglik.true), 2)
  expect_lt(abs(logmeanexp(asifir.loglik) - loglik.true), 2)

})


