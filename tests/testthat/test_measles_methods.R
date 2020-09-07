library(spatPomp)
context("test methods on Measles")
doParallel::registerDoParallel(3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1)
measles3 <- girfd_measles(U=3, N=20, Np = 500, Nguide = 50, lookahead = 1)
girf_loglik <- replicate(10,logLik(girf(measles3)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ABF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1)
measles3 <- abfd_measles(U=3,
                         N=20,
                         Nrep = 100,
                         Np = 50,
                         nbhd = function(object, time, unit) {
                           nbhd_list <- list()
                           if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
                           if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
                           return(nbhd_list)
                         })
abf_loglik <- replicate(10,logLik(abf(measles3)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ABFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1)
measles3 <- abfird_measles(U=3,
                            N=20,
                            Nrep = 100,
                            Np = 50,
                            nbhd = function(object, time, unit) {
                              nbhd_list <- list()
                              if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
                              if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
                              return(nbhd_list)
                            })
abfir_loglik <- replicate(10,logLik(abfir(measles3)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from pfilter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pfilter_loglik <- replicate(10,logLik(pfilter(measles3,
                                            Np = 10000
                                            )))

test_that("ABF, ABFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - logmeanexp(pfilter_loglik)), 15)
  expect_lt(abs(logmeanexp(abf_loglik) - logmeanexp(pfilter_loglik)), 15)
  expect_lt(abs(logmeanexp(abfir_loglik) - logmeanexp(pfilter_loglik)), 15)

})


