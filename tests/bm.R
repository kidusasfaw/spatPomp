png(filename="bm-%02d.png",res=100)
library(spatPomp)


# For CRAN tests, need to limit to two cores
# For covr, needs to be single core (https://github.com/r-lib/covr/issues/227)
doParallel::registerDoParallel(1)
doRNG::registerDoRNG(2)

b_model <- bm(U=2,N=2) 

## ------------------------------------------------------------
## The bm model provides a simple example to test other methods.
## First, we test the filtering methods
## ____________________________________________________________

##
## pfilter tested on bm
##

b_pf <- pfilter(b_model,Np=10)
paste("bm pfilter loglik: ",logLik(b_pf))

##
## abf tested on bm. abf uses parallelization, so we also test that
##

b_bag_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

b_abf <- abf(b_model,Nrep=3,Np=10, nbhd = b_bag_nbhd)
paste("bm abf loglik: ",logLik(b_abf))

##
## abfir tested on bm
##

b_abfir <- abfir(b_model, Nrep = 3, Np = 10, nbhd = b_bag_nbhd)
paste("bm abfir loglik: ",logLik(b_abfir))

##
## bpfilter tested on bm
##

b_bpfilter <- bpfilter(b_model, Np = 10, block_size = 1)
paste("bm bpfilter loglik: ",logLik(b_bpfilter))

##
## enkf tested on bm
##

b_enkf <- enkf(b_model, Np = 10)
paste("bm enkf loglik: ",logLik(b_enkf))

##
## girf tested on bm, both moment and bootstrap methods
##

b_girf_mom <- girf(b_model,Np = 10,lookahead = 1,Nguide = 10,
  kind = 'moment')
paste("bm girf loglik, moment guide: ",logLik(b_girf_mom))

b_girf_boot <- girf(b_model,Np = 10,lookahead = 1,Nguide = 10,
  kind = 'bootstrap')
paste("bm girf loglik, bootstrap guide: ",logLik(b_girf_boot))


## ------------------------------------------------------------
## Now, we test the inference methods
## ____________________________________________________________

b_rw.sd <- rw.sd(rho=0.02,X1_0=ivp(0.02))

##
## igirf on bm
##
## A call to igirf using the moment-based guide function can test compiled
## code for eunit_measure, munit_measure, vunit_measure, dunit_measure,
## runit_measure, rprocess, skeleton, rinit and partrans. 
##
## we test both geometric and hyperbolic cooling

b_igirf_geom <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=10,
  Ninter = 2,
  lookahead = 1,
  Nguide = 5,
  kind = 'moment',
  verbose = FALSE
)
paste("bm igirf loglik, geometric cooling, verbose=F: ",logLik(b_igirf_geom))

b_igirf_hyp <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  Np=10,
  Ninter = 2,
  lookahead = 1,
  Nguide = 5,
  kind = 'moment',
  verbose = TRUE
)
paste("bm igirf loglik, hyperbolic cooling, verbose=T: ",logLik(b_igirf_hyp))

##
## ienkf on bm, with geometric and hyperbolic cooling
##

b_ienkf_geom <- ienkf(b_model,
  Nenkf=2,
  Np = 10,
  rw.sd=b_rw.sd,
  cooling.type="geometric",
  cooling.fraction.50 = 0.5,
  verbose=FALSE
)
paste("bm ienkf loglik, geometric cooling, verbose=F: ",logLik(b_ienkf_geom))

b_ienkf_hyp <- ienkf(b_model,
  Nenkf=2,
  Np = 10,
  rw.sd=b_rw.sd,
  cooling.type="hyperbolic",
  cooling.fraction.50 = 0.5,
  verbose=TRUE
)
paste("bm ienkf loglik, hypoerbolic cooling, verbose=T: ",logLik(b_ienkf_hyp))

##
## iubf on bm, with geometric and hyperbolic cooling
##

b_iubf_geom <- iubf(b_model,
  Nubf = 2,
  Nrep_per_param = 3,
  Nparam = 3,
  nbhd = b_bag_nbhd,
  prop = 0.8,
  rw.sd =b_rw.sd,
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  verbose=FALSE
)
paste("bm iubf loglik, geometric cooling, verbose=F: ",logLik(b_iubf_geom))

b_iubf_hyp <- iubf(b_model,
  Nubf = 2,
  Nrep_per_param = 3,
  Nparam = 3,
  nbhd = b_bag_nbhd,
  prop = 0.8,
  rw.sd =b_rw.sd,
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  verbose=TRUE
)
paste("bm ienkf loglik, hyperbolic cooling, verbose=T: ",logLik(b_iubf_hyp))

## --------------------------------------------
## using bm to test simulate and plot
## ____________________________________________

b_sim1 <- simulate(b_model,nsim=2,format='data.frame')
b_sim2 <- simulate(b_model,nsim=2,format='data.frame',include.data=TRUE)
b_sim3 <- simulate(b_model,nsim=2,format='spatPomps')

plot(b_model,type="l",log=FALSE)
b_sim3v2 <- b_sim3[[1]]
b_sim3v2@data <- exp(b_sim3v2@data)
plot(b_sim3v2,type="l",log=TRUE)
plot(b_sim3[[2]],type="h")

## --------------------------------------------
## using bm to test spatPomp workhorse functions, extending pomp:
## vunit_measure, eunit_measure, munit_measure, dunit_measure
##
## these are tested implicitly in the methods, but here is
## a more direct test
## ____________________________________________


b_s <- states(b_model)[,1,drop=FALSE]
dim(b_s) <- c(dim(b_s),1)
dimnames(b_s) <- list(variable=dimnames(states(b_model))[[1]], rep=NULL)
b_p <- coef(b_model)
dim(b_p) <- c(length(b_p),1)
dimnames(b_p) <- list(param=names(coef(b_model)))

vunit_measure(b_model, x=b_s, unit=2, time=1, params=b_p)

eunit_measure(b_model, x=b_s, unit=2, time=1, params=b_p)

b_array.params <- array(b_p,
  dim = c(length(b_p),length(unit_names(b_model)), 1, 1),
  dimnames = list(params = rownames(b_p)))
b_vc <- c(4, 9) # this should have length equal to the number of units
dim(b_vc) <- c(length(b_vc), 1, 1)

munit_measure(b_model, x=b_s, vc=b_vc, Np=1, unit = 1, time=1,
  params=b_array.params)

dunit_measure(b_model, y=obs(b_model)[,1,drop=FALSE],
  x=b_s, unit=1, time=1, params=b_p)

## --------------------------------------------
## using bm to test edge cases and utility functions
## perhaps only of technical interest
## ____________________________________________

print(b_model)


dev.off()

