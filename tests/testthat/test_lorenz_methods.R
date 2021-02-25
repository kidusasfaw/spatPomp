library(spatPomp)
context("test methods on Lorenz")

# For CRAN tests, need to limit to two cores
# https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions
chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
  # use 2 cores for CRAN
  num_workers <- 2L
} else {
  # use all cores when testing
  num_workers <- parallel::detectCores()
}
num_workers <- 2L
doParallel::registerDoParallel(num_workers)

# create the Lorenz object
set.seed(1)
lorenz4_test <- lorenz(U=4, N=20, delta_t=0.01, delta_obs=1)

# Parameter inference test
start_params <- c('F' = 6, 'sigma' = 0.5, 'tau' = 0.5, "X1_0"=0, "X2_0"=0,
                   "X3_0"=0, "X4_0"=0.01)
## IGIRF
igirf_lookahead <- 1
igirf_ninter <- length(unit_names(lorenz4_test))
igirf_np <- 1000
igirf_nguide <- 40
igirf_ngirf <- 10

igirf_out <- igirf(lorenz4_test,
                   Ngirf = igirf_ngirf,
                   params=start_params,
                   rw.sd = rw.sd(F=0.02, sigma=0.02, tau=0.02,
                                 X1_0=ivp(0),X2_0=ivp(0), X3_0=ivp(0),X4_0=ivp(0)
                   ),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=igirf_np,
                   Ninter = igirf_ninter,
                   lookahead = igirf_lookahead,
                   Nguide = igirf_nguide,
                   kind = 'bootstrap',
                   verbose = FALSE
)

## IUBF
iubf_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(unit>2) nbhd_list <- c(nbhd_list, list(c(unit-2, time)))
  return(nbhd_list)
}
iubf_nubf <- 10
iubf_nrep_per_param <- 1000
iubf_nparam <- 50
iubf_prop <- 0.80

iubf(lorenz4_test,
     Nabf = iubf_nubf,
     Nrep_per_param = iubf_nrep_per_param,
     Nparam = iubf_nparam,
     nbhd = iubf_nbhd,
     params=start_params,
     prop = iubf_prop,
     rw.sd = rw.sd(F=0.02, sigma=0.02, tau=0.02,
                   X1_0=ivp(0),X2_0=ivp(0), X3_0=ivp(0),X4_0=ivp(0)
     ),
     cooling.type = "geometric",
     cooling.fraction.50 = 0.5,
     verbose=FALSE
) -> iubf_out

mif2_out <- mif2(lorenz4_test,
                 params=start_params,
                 Nmif = 20,
                 rw.sd = rw.sd(F=0.02, sigma=0.02, tau=0.02,
                               X1_0=0, X2_0=0, X3_0=0, X4_0=0
                 ),
                 cooling.type = 'geometric',
                 cooling.fraction.50 = 0.5,
                 Np = 1000
)

test_that("IGIRF and IUBF produce estimates that are not far from IF2 for low dimensions", {
  expect_lt(abs(logLik(igirf_out) - logLik(mif2_out)), 20)
  expect_lt(abs(logLik(iubf_out) - logLik(mif2_out)), 35)
})

# Test filtering algorithms
## GIRF
girf_loglik <- replicate(3,logLik(girf(lorenz4_test,
                    Np = 500,
                    lookahead = 1,
                    Nguide = 40
                    )))

## ABF
abf_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(time>2) nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(unit>2) nbhd_list <- c(nbhd_list, list(c(unit-2, time)))

  return(nbhd_list)
}

abf_loglik <- replicate(n=5,
                        expr=logLik(
                          abf(lorenz4_test,
                              Nrep = 500,
                              Np = 20,
                              nbhd = abf_nbhd)))

## ABFIR
abfir_loglik <- replicate(5,
                          logLik(
                            abfir(lorenz4_test,
                                  Nrep = 500,
                                  Np=20,
                                  Nguide = 40,
                                  nbhd = abf_nbhd)))

## PFILTER
pfilter_loglik <- replicate(10,logLik(pfilter(lorenz4_test,
                                            Np = 10000
                                            )))

test_that("ABF, ABFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - logmeanexp(pfilter_loglik)), 10)
  expect_lt(abs(logmeanexp(abf_loglik) - logmeanexp(pfilter_loglik)), 50)
  expect_lt(abs(logmeanexp(abfir_loglik) - logmeanexp(pfilter_loglik)), 30)
})


