png(filename="bm-%02d.png",res=100)
library(spatPomp)
i <- 1
U <- switch(i,2,10)
N <- switch(i,2,10)
Np <- switch(i,10,100)

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

# if(.Platform$OS.type != "windows")
#   doParallel::registerDoParallel(num_workers)

# For covr, needs to be single core (https://github.com/r-lib/covr/issues/227)

# CRAN win-builder test fails in foreach for iubf when using a single
# core registered with 
###  doParallel::registerDoParallel(1)
# so run without registering parallel backend at all
# this generates an R warning
# Warning message:
# executing %dopar% sequentially: no parallel backend registered 
# but that is not a major problem

set.seed(2)
## doRNG::registerDoRNG(2)
## using doRNG with 1 core leads to warnings: it seems to make
## foreach confused about whether it is running in parallel or not.

b_model <- bm(U=U,N=N) 

## ------------------------------------------------------------
## The bm model provides a simple example to test other methods.
## First, we test the filtering methods
## ____________________________________________________________

##
## exact likelihood via the Kalman filter
##

paste("bm kalman filter loglik: ",round(bm_kalman_logLik(b_model),10))

##
## pfilter tested on bm
##

b_pf <- pfilter(b_model,Np=Np)
paste("bm pfilter loglik: ",round(logLik(b_pf),10))

##
## abf tested on bm. abf uses parallelization, so we also test that
##

b_bag_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

set.seed(7)
b_abf <- abf(b_model,Nrep=2,Np=Np, nbhd = b_bag_nbhd)
paste("bm abf loglik: ",round(logLik(b_abf),10))

set.seed(7)
b_abf_repeat <- abf(b_abf)
paste("check abf on abfd_spatPomp: ",
  logLik(b_abf_repeat)==logLik(b_abf))

##
## abfir tested on bm
##

b_abfir <- abfir(b_model, Nrep = 2, Np = Np, nbhd = b_bag_nbhd)
paste("bm abfir loglik: ",round(logLik(b_abfir),10))

##
## bpfilter tested on bm
##

set.seed(5)
b_bpfilter <- bpfilter(b_model, Np = Np, block_size = 1)
paste("bm bpfilter loglik: ",round(logLik(b_bpfilter),10))
set.seed(5)
b_bpfilter_repeat <- bpfilter(b_bpfilter)
paste("check bpfilter on bpfilterd_spatPomp: ",
  logLik(b_bpfilter)==logLik(b_bpfilter_repeat))

##
## enkf tested on bm
##

b_enkf <- enkf(b_model, Np = Np)
paste("bm enkf loglik: ",round(logLik(b_enkf),10))

##
## girf tested on bm, both moment and bootstrap methods
##

b_girf_mom <- girf(b_model,Np = floor(Np/2),lookahead = 1,
  Nguide = floor(Np/2),
  kind = 'moment',Ninter=2)
paste("bm girf loglik, moment guide: ",round(logLik(b_girf_mom),10))

set.seed(0)
b_girf_boot <- girf(b_model,Np = floor(Np/2),lookahead = 1,
  Nguide = floor(Np/2),
  kind = 'bootstrap',Ninter=2)
paste("bm girf loglik, bootstrap guide: ",round(logLik(b_girf_boot),10))

set.seed(0)
b_girf_boot_repeat <- girf(b_girf_boot)
paste("check girf on girfd_spatPomp: ",
  logLik(b_girf_boot)==logLik(b_girf_boot_repeat))

## ------------------------------------------------------------
## Now, we test the inference methods
## ____________________________________________________________

b_rw.sd <- rw_sd(rho=0.02,X1_0=ivp(0.02))

##
## igirf on bm
##
## A call to igirf using the moment-based guide function can test compiled
## code for eunit_measure, munit_measure, vunit_measure, dunit_measure,
## runit_measure, rprocess, skeleton, rinit and partrans. 
##
## we test both geometric and hyperbolic cooling

set.seed(1)
b_igirf_geom <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=Np,
  Ninter = 2,
  lookahead = 1,
  Nguide = 5,
  kind = 'moment',
  verbose = FALSE
)
paste("bm igirf loglik, geometric cooling, verbose=F: ",round(logLik(b_igirf_geom),10))

set.seed(1)
b_igirf_geom_repeat <- igirf(b_igirf_geom,params=coef(b_model))
paste("check igirf on igirfd_spatPomp: ",
  logLik(b_igirf_geom)==logLik(b_igirf_geom_repeat))


b_igirf_hyp <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  Np=floor(Np/2),
  Ninter = 2,
  lookahead = 1,
  Nguide = floor(Np/2),
  kind = 'moment',
  verbose = TRUE
)
paste("bm igirf loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_igirf_hyp),10))

plot(b_igirf_geom) -> b_igirf_plot
head(b_igirf_plot$data)

##
## ienkf on bm, with geometric and hyperbolic cooling
##

b_ienkf_geom <- ienkf(b_model,
  Nenkf=2,
  Np = Np,
  rw.sd=b_rw.sd,
  cooling.type="geometric",
  cooling.fraction.50 = 0.5,
  verbose=FALSE
)
paste("bm ienkf loglik, geometric cooling, verbose=F: ",round(logLik(b_ienkf_geom),10))

b_ienkf_hyp <- ienkf(b_model,
  Nenkf=2,
  Np = Np,
  rw.sd=b_rw.sd,
  cooling.type="hyperbolic",
  cooling.fraction.50 = 0.5,
  verbose=TRUE
)
paste("bm ienkf loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_ienkf_hyp),10))

##
## iubf on bm, with geometric and hyperbolic cooling
##

b_iubf_geom <- iubf(b_model,
  Nubf = 2,
  Nrep_per_param = floor(Np/2),
  Nparam = floor(Np/2),
  nbhd = b_bag_nbhd,
  prop = 0.8,
  rw.sd =b_rw.sd,
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  verbose=FALSE
)
paste("bm iubf loglik, geometric cooling, verbose=F: ",round(logLik(b_iubf_geom),10))

b_iubf_hyp <- iubf(b_model,
  Nubf = 2,
  Nrep_per_param = floor(Np/2),
  Nparam = floor(Np/2),
  nbhd = b_bag_nbhd,
  prop = 0.8,
  rw.sd =b_rw.sd,
#  cooling.type = "hyperbolic",
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  verbose=TRUE
)
paste("bm ienkf loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_iubf_hyp),10))

## --------------------------------------------
## using bm to test simulate and plot
## ____________________________________________

b_sim1 <- simulate(b_model,nsim=2,format='data.frame')
head(b_sim1,10)
b_sim2 <- simulate(b_model,nsim=2,format='data.frame',include.data=TRUE)
head(b_sim2,10)
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
b_y <- obs(b_model)[,1,drop=FALSE]

vunit_measure(b_model, x=b_s, unit=2, time=1, params=b_p)

eunit_measure(b_model, x=b_s, unit=2, time=1, params=b_p)

b_array.params <- array(b_p,
  dim = c(length(b_p),length(unit_names(b_model)), 1, 1),
  dimnames = list(params = rownames(b_p)))
b_vc <- c(4, 9) # this should have length equal to the number of units
dim(b_vc) <- c(length(b_vc), 1, 1)

munit_measure(b_model, x=b_s, vc=b_vc, Np=1, unit = 1, time=1,
  params=b_array.params)

dunit_measure(b_model, y=b_y,
  x=b_s, unit=1, time=1, params=b_p)

runit_measure(b_model, x=b_s, unit=2, time=1, params=b_p)

vec_rmeasure(b_model,x=b_s,time=1, params=b_p)

## --------------------------------------------
## using bm to test edge cases and utility functions
## perhaps only of technical interest
## ____________________________________________

print(b_model)

# check how u is treated by dunit_measure, runit_measure, eunit_measure,
# vunit_measure and munit_measure. this should output unit-1 to
# be consistent with Csnippet indexing.

b_u <- spatPomp(b_model,
  dunit_measure=spatPomp_Csnippet("lik=u;"),
  eunit_measure=spatPomp_Csnippet("ey=u;"),
  munit_measure=spatPomp_Csnippet("M_tau=u;"),
  vunit_measure=spatPomp_Csnippet("vc=u;"),
  runit_measure=spatPomp_Csnippet("Y=u;")  
)

vunit_measure(b_u, x=b_s, unit=2, time=1, params=b_p)
eunit_measure(b_u, x=b_s, unit=2, time=1, params=b_p)
munit_measure(b_u, x=b_s, vc=b_vc, Np=1, unit = 2, time=1,params=b_array.params)
dunit_measure(b_u, y=b_y,x=b_s, unit=2, time=1, params=b_p)
runit_measure(b_u, x=b_s, unit=2, time=1, params=b_p)

dev.off()

## test spatPomp_Csnippet variable construction
spatPomp_Csnippet("lik=u;",unit_statenames="A",unit_obsnames=c("B","C"), unit_covarnames="D",
  unit_ivpnames="E",unit_paramnames="F",unit_vfnames="G")

## --------------------------------------------
## using bm to test spatPomp() replacement functionality
## ____________________________________________

b_rep1 <- spatPomp(b_model,params=coef(b_model))
for(slt in slotNames(b_model)) if(!identical(slot(b_model,slt),slot(b_rep1,slt))) print(slt)

# test an error message
try(spatPomp(data=as.data.frame(b_model),units=NULL),outFile=stdout())

# test parameter replacement
b_rep2 <- spatPomp(b_model,params=coef(b_model)+1)
if(!identical(coef(b_rep2),coef(b_model)+1)) stop('problem with parameter replacement')





