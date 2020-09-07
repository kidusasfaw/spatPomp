library(spatPomp)
context("test methods on Lorenz")
doParallel::registerDoParallel(3)

# create the Lorenz object
set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)
lorenz5_test <- lorenz5
coef(lorenz5_test) <- c('F' = 6, 'sigma' = 0.5, 'tau' = 0.5, "X1_0"=0, "X2_0"=0,
                   "X3_0"=0, "X4_0"=0, "X5_0"=0)
ienkf_Nenkf = 50
ienkf_np = 1000
ienkf_out <- ienkf(lorenz5_test,
                   Nenkf = ienkf_Nenkf,
                   rw.sd = rw.sd(
                     F=1, sigma=1, tau=1, X1_0=0.0, X2_0=0.0,
                     X3_0=0.0,X4_0=0.0,X5_0=0.0),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=ienkf_np)
mif2_out <- mif2(lorenz5_test,
                 Nmif = 100,
                 rw.sd = rw.sd(F=0.02, sigma=0.02, tau=0.02, X1_0=0, X2_0=0,
                               X3_0=0, X4_0=0, X5_0=0),
                 cooling.type = 'geometric',
                 cooling.fraction.50 = 0.5,
                 Np = 1000)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf_loglik <- replicate(10,logLik(girf(lorenz5,
                    Np = 100,
                    lookahead = 1,
                    Nguide = 50
                    )))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ABF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abf_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(unit>2) nbhd_list <- c(nbhd_list, list(c(unit-2, time)))

  return(nbhd_list)
}

abf_loglik <- replicate(n=10,
                        expr=logLik(
                          abf(lorenz5,
                              Nrep = 100,
                              Np = 20,
                              nbhd = abf_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ABFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abfir_loglik <- replicate(10,logLik(abfir(lorenz5,
                        Nrep = 100,
                        Np=20,
                        nbhd = abf_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from pfilter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pfilter_loglik <- replicate(10,logLik(pfilter(lorenz5,
                                            Np = 10000
                                            )))

test_that("ABF, ABFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - logmeanexp(pfilter_loglik)), 10)
  expect_lt(abs(logmeanexp(abf_loglik) - logmeanexp(pfilter_loglik)), 50)
  expect_lt(abs(logmeanexp(abfir_loglik) - logmeanexp(pfilter_loglik)), 30)

})


