library(spatPomp)
context("test senkf2 on Lorenz")
doParallel::registerDoParallel(3)

# create the Lorenz object
set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)

s2 <- senkf2(lorenz5, Np = 1000)

set.seed(1)
lorenz5 <- lorenz(U=5, N=20, dt=0.01, dt_obs=1)
lorenz_h <- function(state.vec, param.vec){
  # extract the indices of the measured components
  ix<-grep('X',names(state.vec))
  # return the measured components
  state.vec[ix]
}

s <- senkf(lorenz5,
      Np = 1000,
      h = lorenz_h,
      R = diag((coef(lorenz5)["tau"])^2,
               nrow = length(spat_units(lorenz5))))
