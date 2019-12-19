library(spatPomp)
context("test methods on Lorenz")
doParallel::registerDoParallel(3)

# create the Lorenz object
set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf_loglik <- replicate(10,logLik(girf(lorenz5,
                    Np = 100,
                    lookahead = 1,
                    Nguide = 50
                    )))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asif_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(unit>2) nbhd_list <- c(nbhd_list, list(c(unit-2, time)))

  return(nbhd_list)
}

asif_loglik <- replicate(10,logLik(asif(lorenz5,
                           islands = 100,
                           Np = 20,
                           nbhd = asif_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asifir_loglik <- replicate(10,logLik(asifir(lorenz5,
                        islands = 100,
                        Np=20,
                        nbhd = asif_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from pfilter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pfilter_loglik <- replicate(10,logLik(pfilter(lorenz5,
                                            Np = 10000
                                            )))

test_that("ASIF, ASIFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - logmeanexp(pfilter_loglik)), 10)
  expect_lt(abs(logmeanexp(asif_loglik) - logmeanexp(pfilter_loglik)), 30)
  expect_lt(abs(logmeanexp(asifir_loglik) - logmeanexp(pfilter_loglik)), 30)

})


