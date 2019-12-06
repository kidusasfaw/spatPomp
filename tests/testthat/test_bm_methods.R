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

fun_to_optim <- function(cf){
  rootQ = cf["rho"]^dmat * cf["sigma"]
  -pomp2:::kalmanFilter(
    t=1:N,
    y=obs(bm3),
    X0=rinit(bm3),
    A=diag(length(spat_units(bm3))),
    Q=rootQ%*%rootQ,
    C=diag(1,nrow=nrow(dmat)),
    R=diag(cf["tau"]^2, nrow=nrow(dmat))
  )$loglik
}
mle <- optim(coef(bm3), fun_to_optim)
kfll_mle <- mle$value
kfll_mle

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   igirf starting from arbitrary parameter set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
igirf.lookahead <- 1
igirf.ninter <- length(spat_units(bm3))
igirf.np <- 800
igirf.nguide <- 40
igirf.ngirf <- 50
coef(bm3) <- c("rho" = 0.7, "sigma"=0.8, "tau"=0.2, "X1_0"=0, "X2_0"=0,
                "X3_0"=0, "X4_0"=0, "X5_0"=0,"X6_0"=0, "X7_0"=0, "X8_0"=0)
igirf.out <- igirf(bm3, Ngirf = igirf.ngirf,
                   rw.sd = rw.sd(rho=0.02, sigma=0.02, tau=0.02, X1_0=0.02,
                                 X2_0=0.02, X3_0=0.02, X4_0=0.02, X5_0=0.02,X6_0=0.02, X7_0=0.02,
                                 X8_0=0.02),
                   cooling.type = "geometric",
                   cooling.fraction.50 = 0.5,
                   Np=igirf.np,
                   Ninter = igirf.ninter,
                   lookahead = igirf.lookahead,
                   Nguide = igirf.nguide)


