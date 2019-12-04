library(spatPomp)
context("test methods on Lorenz")
doParallel::registerDoParallel(3)

# create the Lorenz object
lorenz5 <- lorenz(U=5, N=100, dt=0.01, dt_obs=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf.loglik <- replicate(10,logLik(girf(lorenz5,
                    Np = 100,
                    Ninter = length(spat_units(lorenz5)),
                    lookahead = 1,
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

asif.loglik <- replicate(10,logLik(asif(lorenz5,
                           islands = 50,
                           Np = 20,
                           nbhd = asif.nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asifir.loglik <- replicate(10,logLik(asifir(lorenz5,
                        islands = 50,
                        Np=20,
                        nbhd = asif.nbhd,
                        Ninter = length(spat_units(lorenz5)))))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from pfilter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pfilter.loglik <- replicate(10,logLik(pfilter(lorenz5,
                                            Np = 2000
                                            )))

test_that("ASIF, ASIFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf.loglik) - logmeanexp(pfilter.loglik)), 2)
  expect_lt(abs(logmeanexp(asif.loglik) - logmeanexp(pfilter.loglik)), 2)
  expect_lt(abs(logmeanexp(asifir.loglik) - logmeanexp(pfilter.loglik)), 2)

})


