library(spatPomp)
context("test methods on simple Brownian motion")

doParallel::registerDoParallel(3)
# create the BM object
set.seed(1)
U = 8; N = 10
bm8 <- bm(U = U, N = N)
bpfilter_loglik <- replicate(10, bpfilter(bm8, Np = 1000, num_partitions = 4)@loglik)
pfilter_loglik <- replicate(10, pfilter(bm8, Np = 100)@loglik)
# compute distance matrix to compute true log-likelihood
dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
dmat <- matrix(0,U,U)
for(u in 1:U) {
  for(v in 1:U) {
    dmat[u,v] <- dist(u,v)
  }
}

# compute the true log-likelihood
rootQ = coef(bm8)["rho"]^dmat * coef(bm8)["sigma"]
loglik_true <- pomp:::kalmanFilter(
  t=1:N,
  y=obs(bm8),
  X0=rinit(bm8),
  A= diag(length(spat_units(bm8))),
  Q=rootQ%*%rootQ,
  C=diag(1,nrow=nrow(dmat)),
  R=diag(coef(bm8)["tau"]^2, nrow=nrow(dmat))
)$loglik

fun_to_optim <- function(cf){
  rootQ = cf["rho"]^dmat * cf["sigma"]
  -pomp2:::kalmanFilter(
    t=1:N,
    y=obs(bm8),
    X0=rinit(bm8),
    A=diag(length(spat_units(bm8))),
    Q=rootQ%*%rootQ,
    C=diag(1,nrow=nrow(dmat)),
    R=diag(cf["tau"]^2, nrow=nrow(dmat))
  )$loglik
}
mle <- optim(coef(bm8), fun_to_optim)
kfll_mle <- mle$value
kfll_mle

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf_loglik <- replicate(10,logLik(girf(bm8,
                    Np = 500,
                    lookahead = 1,
                    Nguide = 50
                    )))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from GIRF with lookahead > 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girf_loglik_l2 <- replicate(10,logLik(girf(bm3,
                                        Np = 500,
                                        lookahead = 2,
                                        Nguide = 50
)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asif_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

asif_loglik <- replicate(10,logLik(asif(bm8,
                           islands = 100,
                           Np = 50,
                           nbhd = asif_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from ASIFIR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
asifir_loglik <- replicate(10,logLik(asifir(bm3,
                        islands = 100,
                        Np=50,
                        nbhd = asif_nbhd)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from EnKF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
genkf_loglik <- replicate(10,logLik(genkf(bm8, Np = 1000)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   log-likelihood estimate from bpfilter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



test_that("ASIF, ASIFIR, GIRF all yield close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(asif_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(asifir_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(genkf_loglik) - loglik_true), 3)
  expect_lt(abs(logmeanexp(asifir_loglik) - loglik_true), 3)


})

test_that("GIRF with lookahead >= 2 yields close to true log-likelihood estimates", {
  expect_lt(abs(logmeanexp(girf_loglik_l2) - loglik_true), 3)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   igirf starting from arbitrary parameter set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
igirf_lookahead <- 1
igirf_ninter <- length(spat_units(bm3))
igirf_np <- 800
igirf_nguide <- 40
igirf_ngirf <- 30
coef(bm3) <- c("rho" = 0.7, "sigma"=0.8, "tau"=0.2, "X1_0"=0, "X2_0"=0,
                "X3_0"=0, "X4_0"=0, "X5_0"=0,"X6_0"=0, "X7_0"=0, "X8_0"=0)
igirf_out1 <- igirf(bm3, Ngirf = igirf_ngirf,
                   rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=0.02,
                                 X2_0=0.02, X3_0=0.02, X4_0=0.02, X5_0=0.02,X6_0=0.02, X7_0=0.02,
                                 X8_0=0.02),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=igirf_np,
                   Ninter = igirf_ninter,
                   lookahead = igirf_lookahead,
                   Nguide = igirf_nguide)
igirf_out2 <- igirf(bm3, Ngirf = igirf_ngirf,
                   rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=0.02,
                                 X2_0=0.02, X3_0=0.02, X4_0=0.02, X5_0=0.02,X6_0=0.02, X7_0=0.02,
                                 X8_0=0.02),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=igirf_np,
                   Ninter = igirf_ninter,
                   lookahead = igirf_lookahead,
                   Nguide = igirf_nguide,
                   method = 'adams')


test_that("IGIRF produces estimates that are not far from the MLE", {
  expect_lt(abs(logLik(igirf_out1) - (-kfll_mle)), 20)
  expect_lt(abs(logLik(igirf_out2) - (-kfll_mle)), 20)

})

