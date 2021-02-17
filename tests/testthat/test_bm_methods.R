library(spatPomp)
context("Test methods on simple Brownian motion")

doParallel::registerDoParallel(3)
# create the BM object
set.seed(2)
U = 10; N = 5
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
                  "X1_0"=0, "X2_0"=0, "X3_0"=0, "X4_0"=0, "X5_0"=0,
                  "X6_0"=0, "X7_0"=0, "X8_0"=0, "X9_0"=0, "X10_0"=0)

## IEnKF
ienkf_np <- 1000
ienkf_Nenkf <- 20
ienkf_out <- ienkf(bm_obj,
                   Nenkf = ienkf_Nenkf,
                   params = start_params,
                   rw.sd = rw.sd(
                     rho=0.02, sigma=0.02, tau=0.02,
                     X1_0=0.0, X2_0=0.0, X3_0=0.0, X4_0=0.0, X5_0=0.0,
                     X6_0=0.0, X7_0=0.0, X8_0=0.0, X9_0=0.0, X10_0=0.0),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=ienkf_np)

## IGIRF
igirf_lookahead <- 1
igirf_ninter <- length(unit_names(bm_obj))
igirf_np <- 800
igirf_nguide <- 40
igirf_ngirf <- 30
### Use moment-matching approach
igirf_out1 <- igirf(bm_obj, Ngirf = igirf_ngirf,
                    params=start_params,
                    rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=0.02,
                                  X2_0=0.02, X3_0=0.02, X4_0=0.02, X5_0=0.02,X6_0=0.02, X7_0=0.02,
                                  X8_0=0.02, X9_0=0.02, X10_0=0.02),
                    cooling.type = "geometric",
                    cooling.fraction.50 = 0.5,
                    Np=igirf_np,
                    Ninter = igirf_ninter,
                    lookahead = igirf_lookahead,
                    Nguide = igirf_nguide,
                    kind = 'moment')
### Use quantile-based approach
igirf_out2 <- igirf(bm_obj, Ngirf = igirf_ngirf,
                    params=start_params,
                    rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=ivp(0.02),
                                  X2_0=ivp(0.02), X3_0=ivp(0.02), X4_0=ivp(0.02), X5_0=ivp(0.02),X6_0=ivp(0.02), X7_0=ivp(0.02),
                                  X8_0=ivp(0.02), X9_0=ivp(0.02), X10_0=ivp(0.02)),
                    cooling.type = "geometric",
                    cooling.fraction.50 = 0.5,
                    Np=igirf_np,
                    Ninter = igirf_ninter,
                    lookahead = igirf_lookahead,
                    Nguide = igirf_nguide,
                    kind = 'bootstrap')

test_that("IGIRF produces estimates that are not far from the MLE", {
  expect_lt(abs(logLik(igirf_out1) - (-kfll_mle)), 20)
  expect_lt(abs(logLik(igirf_out2) - (-kfll_mle)), 20)
})
### IEnKF will struggle with the bm problem
test_that("IEnKF produces rho value that is not far from the MLE", {
  expect_lt(abs(coef(ienkf_out)['rho'] - (mle$par['rho'])), 0.3)
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
girf_loglik <- replicate(n = 10,
                         expr = logLik(
                           girf(bm_obj,
                                Np = 500,
                                lookahead = 1,
                                Nguide = 50
                            )
                          )
)

## GIRF with lookahead > 1
girf_loglik2 <- replicate(n = 5,
                         expr = logLik(
                           girf(bm_obj,
                                Np = 500,
                                lookahead = 2,
                                Nguide = 50
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

abf_loglik <- replicate(n=10,
                         expr = logLik(abf(bm_obj,
                                           Nrep = 800,
                                           Np = 50,
                                           nbhd = abf_nbhd)))

## ABFIR
abfir_loglik <- replicate(n=10,
                          expr = logLik(abfir(bm_obj,
                                                Nrep = 800,
                                                Np=50,
                                                nbhd = abf_nbhd)))

## BPF
bpfilter_loglik <- replicate(10,logLik(bpfilter(bm_obj, Np = 100, block_size = 3)))

test_that("ABF, ABFIR, GIRF, EnKF, BPF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(abf_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(abfir_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(enkf_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(bpfilter_loglik) - loglik_true), 3)
})

test_that("GIRF with lookahead >= 2 yields close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik_l2) - loglik_true), 3)
})
