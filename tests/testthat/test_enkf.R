library(spatPomp)
context("test enkf on Lorenz")
doParallel::registerDoParallel(3)
set.seed(1)
U = 5; N = 50

# create the Lorenz object
set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)
# output from enkf
gl <- enkf(lorenz5, Np = 1000)

# recreate the same Lorenz object
set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)
lorenz_h <- function(state.vec, param.vec){
  # extract the indices of the measured components
  ix<-grep('X',names(state.vec))
  # return the measured components
  state.vec[ix]
}
# output from pomp::enkf
el <- pomp::enkf(lorenz5,
      Np = 1000,
      h = lorenz_h,
      R = diag((coef(lorenz5)["tau"])^2,
               nrow = length(unit_names(lorenz5))))

test_that("pomp::enkf and spatPomp::enkf yield equal log likelihoods when vmeasure is independent of X", {
  expect_equal(logLik(el), logLik(gl))
})
