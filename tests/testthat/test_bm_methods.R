library(spatPomp)
context("test methods on simple Brownian motion")

doParallel::registerDoParallel(3)
# create the BM object
U = 8; N = 10
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

girf.loglik <- replicate(10,logLik(girf(bm3,
                    Np = 500,
                    Ninter = length(spat_units(bm3)),
                    lookahead = 1,
                    Nguide = 50
                    )))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF with lookahead > 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf.loglik.l2 <- replicate(10,logLik(girf(bm3,
                                        Np = 500,
                                        Ninter = length(spat_units(bm3)),
                                        lookahead = 2,
                                        Nguide = 50
)))

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
                           islands = 100,
                           Np = 50,
                           nbhd = asif.nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asifir.loglik <- replicate(10,logLik(asifir(bm3,
                        islands = 100,
                        Np=50,
                        nbhd = asif.nbhd,
                        Ninter = length(spat_units(bm3)))))

test_that("ASIF, ASIFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf.loglik) - loglik.true), 3)
  expect_lt(abs(logmeanexp(asif.loglik) - loglik.true), 3)
  expect_lt(abs(logmeanexp(asifir.loglik) - loglik.true), 3)

})

test_that("GIRF with lookahead >= 2 yields close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf.loglik.l2) - loglik.true), 2)
})


