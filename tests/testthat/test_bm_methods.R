library(spatPomp)
context("Test methods on simple Brownian motion")

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
if(.Platform$OS.type != "windows")
  doParallel::registerDoParallel(num_workers)

# create the BM object
set.seed(2)
U = 2; N = 10
bm_obj <- bm(U = U, N = N)


# compute distance matrix to compute true log-likelihood
dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
dmat <- matrix(0,U,U)
for(u in 1:U) {
  for(v in 1:U) {
    dmat[u,v] <- dist(u,v)
  }
}

# compute the true log-likelihood
rootQ = coef(bm_obj)["rho"]^dmat * coef(bm_obj)["sigma"]
loglik_true <- pomp:::kalmanFilter(
  t=1:N,
  y=obs(bm_obj),
  X0=rinit(bm_obj),
  A= diag(length(unit_names(bm_obj))),
  Q=rootQ%*%rootQ,
  C=diag(1,nrow=nrow(dmat)),
  R=diag(coef(bm_obj)["tau"]^2, nrow=nrow(dmat))
)$loglik

fun_to_optim <- function(cf){
  rootQ = cf["rho"]^dmat * cf["sigma"]
  -pomp:::kalmanFilter(
    t=1:N,
    y=obs(bm_obj),
    X0=rinit(bm_obj),
    A=diag(length(unit_names(bm_obj))),
    Q=rootQ%*%rootQ,
    C=diag(1,nrow=nrow(dmat)),
    R=diag(cf["tau"]^2, nrow=nrow(dmat))
  )$loglik
}
mle <- optim(coef(bm_obj), fun_to_optim)
kfll_mle <- -mle$value

# Test inference algorithms
start_params <- c("rho" = 0.7, "sigma"=0.5, "tau"=0.5,
                  "X1_0"=0, "X2_0"=0, "X3_0"=0
                 )

## IGIRF
igirf_lookahead <- 1
igirf_ninter <- length(unit_names(bm_obj))
igirf_np <- 500
igirf_nguide <- 40
igirf_ngirf <- 10
igirf_out <- igirf(bm_obj, Ngirf = igirf_ngirf,
                    params=start_params,
                    rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02,
                                  X1_0=ivp(0), X2_0=ivp(0), X3_0=ivp(0)
                                 ),
                    cooling.type = "geometric",
                    cooling.fraction.50 = 0.5,
                    Np=igirf_np,
                    Ninter = igirf_ninter,
                    lookahead = igirf_lookahead,
                    Nguide = igirf_nguide,
                    kind = 'moment',
                    verbose = FALSE)

test_that("IGIRF produces estimates that are not far from the MLE", {
  expect_lt(abs(logLik(igirf_out) - kfll_mle), 20)
})

# IUBF
iubf_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}
iubf_nubf <- 10
iubf_nrep_per_param <- 500
iubf_nparam <- 50
iubf_prop <- 0.80

iubf(bm_obj,
     Nubf = iubf_nubf,
     Nrep_per_param = iubf_nrep_per_param,
     Nparam = iubf_nparam,
     nbhd = iubf_nbhd,
     params=start_params,
     prop = iubf_prop,
     rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=ivp(0),
                   X2_0=ivp(0), X3_0=ivp(0)
     ),
     cooling.type = "geometric",
     cooling.fraction.50 = 0.5,
     verbose=FALSE
) -> iubf_out

test_that("IUBF produces estimates that are not far from the MLE", {
  expect_lt(abs(logLik(iubf_out) - kfll_mle), 20)
})

# Test filtering algorithms

## EnKF
enkf_loglik <- replicate(n = 10,
                         expr = logLik(
                           enkf(bm_obj,
                                Np = 500
                           )
                         )
)

## GIRF
girf_loglik <- replicate(n = 3,
                         expr = logLik(
                           girf(bm_obj,
                                Np = 500,
                                lookahead = 1,
                                Nguide = 50,
                                kind = 'moment'
                            )
                          )
)

## ABF
abf_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

abf_runs <- 2
abf_loglik <- vector(length = abf_runs)
for(i in seq_len(abf_runs)){
  abf_loglik[i] <- logLik(abf(bm_obj,
                              Nrep = 200,
                              Np = 50,
                              nbhd = abf_nbhd))
}

## ABFIR
abfir_loglik <- vector(length = abf_runs)
for(i in seq_len(abf_runs)){
  abfir_loglik[i] <- logLik(abfir(bm_obj,
                                  Nrep = 200,
                                  Np = 50,
                                  nbhd = abf_nbhd))
}

## BPF
bpfilter_loglik <- replicate(10,logLik(bpfilter(bm_obj, Np = 500, block_size = 1)))

test_that("ABF, ABFIR, GIRF, EnKF, BPF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - loglik_true), 10)
  expect_lt(abs(logmeanexp(abf_loglik) - loglik_true), 10)
  expect_lt(abs(logmeanexp(abfir_loglik) - loglik_true), 10)
  expect_lt(abs(logmeanexp(enkf_loglik) - loglik_true), 10)
  expect_lt(abs(logmeanexp(bpfilter_loglik) - loglik_true), 10)
})
