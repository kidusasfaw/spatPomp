library(spatPomp)
i <- 1
U <- switch(i,2,10)
N <- switch(i,2,10)
Np <- switch(i,5,100)

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
b_model_no_params <-  spatPomp(b_model,params=NULL)
b_model_no_rprocess <- spatPomp(b_model,rprocess=NULL)
b_model_no_eunit_measure <- spatPomp(b_model,eunit_measure=NULL)
b_model_no_vunit_measure <- spatPomp(b_model,vunit_measure=NULL)
b_model_no_runit_measure <- spatPomp(b_model,runit_measure=NULL)
b_model_no_dunit_measure <- spatPomp(b_model,dunit_measure=NULL)
b_model_no_munit_measure <- spatPomp(b_model,munit_measure=NULL)
b_model_t0_equal_t1 <- spatPomp(b_model,t0=1)
b_model5 <- bm(U=U,N=5) 
b_model_with_accumvars <- b_model
b_model_with_accumvars@accumvars <- rownames(states(b_model))

b_model_zero_dmeasure <- spatPomp(b_model,
  dmeasure = spatPomp_Csnippet(
    method="dmeasure",
    unit_statenames="X",
    unit_obsnames="Y",
    code = "
      lik = give_log ? log(0) : 0;
    "
  ),
  dunit_measure = spatPomp_Csnippet("
    lik = give_log ? log(0) : 0;
  ")
)

b_bag_nbhd <- function(object, time, unit) {
  nbhd_list <- list()
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  return(nbhd_list)
}

b_bag_nbhd_lookback1 <- function(object, time, unit) {
  nbhd_list <- list()
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(time>1)  nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  return(nbhd_list)
}

b_bag_nbhd_lookback2 <- function(object, time, unit) {
  nbhd_list <- list()
  if(unit>1) nbhd_list <- c(nbhd_list, list(c(unit-1, time)))
  if(time>1) nbhd_list <- c(nbhd_list, list(c(unit, time-1)))
  if(time>2)  nbhd_list <- c(nbhd_list, list(c(unit, time-2)))
  return(nbhd_list)
}

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

## ---------------------------------------------------------------
## abf tested on bm. abf uses parallelization, so we also test that
## _______________________________________________________________


set.seed(7)
b_abf <- abf(b_model,Nrep=2,Np=Np, nbhd = b_bag_nbhd)
paste("bm abf loglik: ",round(logLik(b_abf),10))

set.seed(7)
b_abf_repeat <- abf(b_abf)
paste("check abf on abfd_spatPomp: ",
  logLik(b_abf_repeat)==logLik(b_abf))

try(abf(b_model))
try(abf(b_model,Nrep=2))
try(abf(b_model,Nrep=2,Np=3))
try(abf(b_model_no_params,Nrep=2,Np=3))
try(abf(b_model,Nrep=2,Np="JUNK"))
try(abf(b_model,Nrep=2,Np=function(n)"JUNK"))
try(abf(b_model,Nrep=2,Np=-1))
try(abf(b_model,Nrep=2,Np=1:42))
param_matrix <- cbind(coef(b_model),coef(b_model))
rownames(param_matrix) <- names(coef(b_model))
try(abf(b_model,Nrep=2,Np=2,params=param_matrix))
abf(b_model_zero_dmeasure,Nrep=2,Np=Np,verbose=TRUE)


## ---------------------------------------------------------------
## abfir tested on bm
## _______________________________________________________________

b_abfir <- abfir(b_model, Nrep = 2, Np = Np, nbhd = b_bag_nbhd)
paste("bm abfir loglik: ",round(logLik(b_abfir),10))

capture.output(
  abfir(b_abfir,verbose=TRUE,accumvars="X1")
) -> b_abfir_out
abfir(b_model, Nrep = 2, Np = Np)
try(abfir(b_model))
try(abfir(b_model,Nrep=2))
try(abfir(b_model,Nrep=2,Np=function(n)"JUNK"))
try(abfir(b_model,Nrep=2,Np=1:10))
try(abfir(b_model,Nrep=2,Np=Np,params=unname(coef(b_model))))
try(abfir(b_model_zero_dmeasure,Nrep = 3, Np = Np))

# test abfir when all particles fail...
# make this happen by setting a high tol
abfir(b_abfir,tol=1000)

## --------------------------------------------------------------
## bpfilter tested on bm
## ______________________________________________________________

set.seed(5)
b_bpfilter <- bpfilter(b_model, Np = Np, block_size = 1)
paste("bm bpfilter loglik: ",round(logLik(b_bpfilter),10))
set.seed(5)
b_bpfilter_repeat <- bpfilter(b_bpfilter)
paste("check bpfilter on bpfilterd_spatPomp: ",
  logLik(b_bpfilter)==logLik(b_bpfilter_repeat))

bpfilter(b_model_t0_equal_t1,Np = Np, block_size = 1,filter_traj=TRUE)

set.seed(5)
b_bpfilter_filter_traj <- bpfilter(b_bpfilter,filter_traj=TRUE)
paste("bpfilter filter trajectory final particle: ")
round(b_bpfilter_filter_traj@filter.traj[,1,],3)

set.seed(5)
b_bpfilter_save_states <- bpfilter(b_bpfilter,save_states=TRUE)
paste("bpfilter final particles: ")
round(b_bpfilter_save_states@saved.states[[N]],3)

set.seed(5)
b_bpfilter_inf <- bpfilter(b_model_zero_dmeasure, Np = Np, block_size = 1)
paste("bm bpfilter loglik, zero measurement: ",
  round(logLik(b_bpfilter_inf),10))

## test bpfilter error messages
try(bpfilter())
try(bpfilter("JUNK"))
try(bpfilter(b_model))
try(bpfilter(b_model,block_list=block_list,block_size=23))
try(bpfilter(b_model,block_list=block_list))
try(bpfilter(b_model,Np=10,block_size=1000))
try(bpfilter(b_bpfilter,block_list=block_list,block_size=23))
try(bpfilter(b_bpfilter,Np=10,block_size=1000))
try(bpfilter(b_bpfilter,Np=-1))
try(bpfilter(b_bpfilter,Np=1:1000))
try(bpfilter(b_bpfilter,Np="JUNK"))
test_params_matrix <- cbind(coef(b_model),coef(b_model),coef(b_model))
try(bpfilter(b_bpfilter,params=test_params_matrix))


## -----------------------------------------------------------------
## enkf tested on bm
## ________________________________________________________________

## test error messages
try(enkf())
try(enkf("JUNK"))
try(enkf(b_model_no_rprocess))
try(enkf(b_model_no_eunit_measure))
try(enkf(b_model_no_vunit_measure))

set.seed(5)
b_enkf <- enkf(b_model, Np = Np)
paste("bm enkf loglik: ",round(logLik(b_enkf),10))


## -----------------------------------------------------------------
## girf on bm: moment and bootstrap methods, followed by error tests
## ________________________________________________________________

set.seed(0)
b_girf_mom <- girf(b_model,Np = floor(Np/2),lookahead = 1,
  Nguide = floor(Np/2),
  kind = 'moment',Ninter=2)
paste("bm girf loglik, moment guide: ",round(logLik(b_girf_mom),10))

## for boostrap girf, we do not set Ninter, to test the default which is Ninter=U
set.seed(0)
b_girf_boot <- girf(b_model,Np = floor(Np/2),lookahead = 1,
  Nguide = floor(Np/2),
  kind = 'bootstrap')
paste("bm girf loglik, bootstrap guide: ",round(logLik(b_girf_boot),10))

set.seed(0)
b_girf_boot_repeat <- girf(b_girf_boot)
paste("check girf on girfd_spatPomp: ",
  logLik(b_girf_boot)==logLik(b_girf_boot_repeat))

## check girf for zero measurement density situations
b_girf_mom_inf <- girf(b_model_zero_dmeasure,Np = floor(Np/2),lookahead = 1,
  Nguide = 3,
  kind = 'moment',Ninter=2)
paste("bm moment girf loglik, zero measurement: ",
  round(logLik(b_girf_mom_inf),10))

set.seed(0)
b_girf_boot_inf <- girf(b_model_zero_dmeasure,Np = floor(Np/2),lookahead = 1,
  Nguide = 3,
  kind = 'bootstrap',Ninter=2)
paste("bm bootstrap girf loglik, zero measurement: ",
  round(logLik(b_girf_boot_inf),10))

print("The following deliver an error message, to test it")
try(girf())
try(girf("JUNK"))
try(girf(b_girf_boot,Np=c(Inf)))
try(girf(b_girf_boot,Np=seq(from=10,length=N+1,by=2)))
try(girf(b_model_no_eunit_measure,kind='moment'))
try(girf(b_model_no_vunit_measure,kind='moment'))
try(girf(b_model_no_rprocess,kind='moment'))
try(girf(b_model,kind='moment'))
try(girf(b_model,kind='moment',Np=5))
try(girf(b_model,kind='moment',Np=5,Nguide=3,tol=1:1000))
try(girf(b_model_no_rprocess,kind='boot'))
try(girf(b_model,kind='boot'))
try(girf(b_model,kind='boot',Np=5))
try(girf(b_model,kind='boot',Np=5,Nguide=3,tol=1:1000))

try(girf(b_model_no_params,Np = 3,lookahead = 1, Nguide = 3,
  kind = 'moment',Ninter=2))
try(girf(b_model_no_params,Np = 3,lookahead = 1, Nguide = 3,
  kind = 'boot',Ninter=2))
  
try(girf(b_model,Np = 1:10,lookahead = 1, Nguide = 3,
  kind = 'moment',Ninter=2))
try(girf(b_model,Np = "JUNK",lookahead = 1, Nguide = 3,
  kind = 'moment',Ninter=2))

girf(b_model_with_accumvars,Np = 3,lookahead = 2, Nguide = 3,
  kind = 'moment',Ninter=2)
girf(b_model_with_accumvars,Np = 3,lookahead = 2, Nguide = 3,
  kind = 'boot',Ninter=2)

# test girf when all particles fail...
# make this happen by setting a high tol
girf(b_girf_mom,tol=1000)



## ------------------------------------------------------------
## Now, we test the inference methods
## ____________________________________________________________

b_rw.sd <- rw_sd(rho=0.02,X1_0=ivp(0.02))

##############################################################
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
  lookahead = 2,
  Nguide = 3,
  kind = 'moment',
  verbose = FALSE
)
paste("bm igirf loglik, geometric cooling, verbose=F: ",round(logLik(b_igirf_geom),5))

set.seed(1)
b_igirf_geom_repeat <- igirf(b_igirf_geom,params=coef(b_model))
paste("check igirf on igirfd_spatPomp: ",
  logLik(b_igirf_geom)==logLik(b_igirf_geom_repeat))

b_igirf_hyp <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  Np=3,
  Ninter = 2,
  lookahead = 2,
  Nguide = floor(Np/2),
  kind = 'moment',
  verbose = TRUE
)
paste("bm igirf loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_igirf_hyp),10))

plot(b_igirf_geom) -> b_igirf_plot
head(b_igirf_plot$data)

set.seed(1)
b_igirf_boot_geom <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  Np=Np,
  Ninter = 2,
  lookahead = 1,
  Nguide = 5,
  kind = 'bootstrap',
  verbose = FALSE
)
paste("bm igirf boot loglik, geometric cooling, verbose=F: ",round(logLik(b_igirf_boot_geom),10))

b_igirf_boot_hyp <- igirf(b_model,
  Ngirf = 2,
  rw.sd = b_rw.sd,
  cooling.type = "hyperbolic",
  cooling.fraction.50 = 0.5,
  Np=3,
  Ninter = 2,
  lookahead = 2,
  Nguide = 3,
  kind = 'bootstrap',
  verbose = TRUE
)
paste("bm igirf boot loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_igirf_hyp),10))


print("The following deliver an error message, to test it")
try(igirf())
try(igirf(data="JUNK"))
try(igirf(b_igirf_boot_geom,Np=c(Inf)))
igirf(b_igirf_boot_geom,Np=3)
try(igirf(b_igirf_boot_geom,Np=NULL))
try(igirf(b_model_no_eunit_measure,kind='moment', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=floor(Np/2), Ninter = 2))
try(igirf(b_model_no_vunit_measure,kind='moment', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=floor(Np/2), Ninter = 2))
try(igirf(b_model_no_rprocess,kind='moment', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=floor(Np/2), Ninter = 2))

try(igirf(b_model))
try(igirf(b_model,Ngirf=2))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=0.5))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=0.5,
  cooling.type="geometric"))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=0.5,
  cooling.type="geometric",Np=4))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=0.5,
  cooling.type="geometric",Np="JUNK",Nguide=4))
try(igirf(b_model,Ngirf="JUNK",rw.sd=b_rw.sd,cooling.fraction.50=0.5,
  cooling.type="geometric",Np=3,Nguide=4))
try(igirf(b_model,Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=1000,
  cooling.type="geometric",Np=3,Nguide=4))

try(igirf(b_model,kind="moment",Ngirf=2,rw.sd=b_rw.sd,cooling.fraction.50=0.5,cooling.type="geometric",Np=3,Nguide=4,tol=-1))


igirf(b_igirf_boot_geom,Np=3,
  .paramMatrix=cbind(coef(b_model),coef(b_model),coef(b_model)))
try(igirf(b_igirf_boot_geom,Np=function(x) 4))
try(igirf(b_igirf_boot_geom,Np=function(x) "JUNK"))
try(igirf(b_igirf_boot_geom,Np=5:15))
try(igirf(b_igirf_boot_geom,tol=-1))


igirf(b_model_with_accumvars,kind='moment', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=3, Ninter = 1)
igirf(b_model_with_accumvars,kind='boot', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=3, Ninter = 1)

igirf(b_model_zero_dmeasure,kind='moment', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=3, Ninter = 1)
igirf(b_model_zero_dmeasure,kind='boot', Ngirf = 2, Nguide=2,
  rw.sd = b_rw.sd, cooling.type = "hyperbolic", cooling.fraction.50 = 0.5,
  Np=3, Ninter = 1)


## ----------------------------------------------------------
## ienkf on bm, with geometric and hyperbolic cooling
## __________________________________________________________

set.seed(55)

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

## test error messages for ienkf
try(ienkf(b_model_no_rprocess))
try(ienkf(b_model, Nenkf="JUNK"))
ienkf(b_model,Nenkf=2,Np = 3,rw.sd=b_rw.sd,cooling.type="geometric",
  cooling.fraction.50 = 0.5,
  .paramMatrix=cbind(coef(b_model),coef(b_model),coef(b_model)))
try(ienkf(b_model,Nenkf=2))
try(ienkf(b_model,Nenkf=2,Np=NULL))
try(ienkf(b_model,Nenkf=2,Np="JUNK"))
try(ienkf(b_model,Nenkf=2,Np = 3))
try(ienkf(b_model,Nenkf=2,Np = 3,rw.sd=b_rw.sd))
try(ienkf(b_model,Nenkf=2,Np = 3,rw.sd=b_rw.sd,cooling.fraction.50 = 1000))
try(ienkf(b_model,Nenkf=2,Np = 3,rw.sd=b_rw.sd,cooling.fraction.50 = 0.5,.indices=1:1000))


## ---------------------------------------------------------------
## iubf on bm, with geometric and hyperbolic cooling
## _______________________________________________________________

set.seed(8)

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
paste("bm iubf loglik, geometric cooling, verbose=F: ",round(logLik(b_iubf_geom),10))

b_iubf_hyp <- iubf(b_model,
  Nubf = 2,
  Nrep_per_param = 3,
  Nparam = 3,
  nbhd = b_bag_nbhd,
  prop = 0.8,
  rw.sd =b_rw.sd,
#  cooling.type = "hyperbolic",
  cooling.type = "geometric",
  cooling.fraction.50 = 0.5,
  verbose=TRUE
)
paste("bm iubf loglik, hyperbolic cooling, verbose=T: ",round(logLik(b_iubf_hyp),10))

try(iubf(b_model))
try(iubf(b_model,Nubf=-1))
try(iubf(b_model,Nubf=2))
try(iubf(b_model,Nubf=2,Nrep_per_param=3))
try(iubf(b_model,Nubf=2,Nrep_per_param=3,rw.sd=b_rw.sd))
try(iubf(b_model,Nubf=2,Nrep_per_param=3,rw.sd=b_rw.sd,cooling.fraction.50=1000))
try(iubf(b_model,Nubf=2,Nrep_per_param=3,rw.sd=b_rw.sd,cooling.fraction.50=0.5))
try(iubf(b_model,Nubf=2,Nrep_per_param=3,rw.sd=b_rw.sd,cooling.fraction.50=0.5,nbhd=b_bag_nbhd))
try(iubf(b_model,Nubf=2,Nrep_per_param=3,rw.sd=b_rw.sd,cooling.fraction.50=0.5,nbhd=b_bag_nbhd,Nparam=3))
try(iubf(b_model,Nubf=2,Nrep_per_param=1,rw.sd=b_rw.sd,cooling.fraction.50=0.5,Nparam=3,nbhd=b_bag_nbhd))
try(iubf(b_model,Nubf=2,Nrep_per_param=1,rw.sd=b_rw.sd,cooling.fraction.50=1000,Nparam=3,nbhd=b_bag_nbhd,prop=0.8))

## max_lookback is not triggered for b_model with N=2
iubf(b_model5,Nubf=2, Nparam = 3,Nrep_per_param=3,nbhd=b_bag_nbhd,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.8)

set.seed(98)

## trigger special case when Nrep_per_param=1
b_iubf_npp1 <- iubf(b_model,Nubf=2, Nparam = 3,Nrep_per_param=1,nbhd=b_bag_nbhd,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.8)
paste("bm iubf loglik with Nrep_per_param = 1 : ",round(logLik(b_iubf_npp1),10))

## trigger special cases when length(def_resample)==0
b_iubf_np1 <- iubf(b_model,Nubf=2, Nparam = 3, Nrep_per_param=1,nbhd=b_bag_nbhd,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0)
paste("bm iubf loglik with length(def_resample)==0: ",round(logLik(b_iubf_np1),10))

## trigger special cases when length(def_resample)==Nparam*prop=1
b_iubf_dr1 <- iubf(b_model,Nubf=2, Nparam = 3,Nrep_per_param=1,nbhd=b_bag_nbhd,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.25)
paste("bm iubf loglik with length(def_resample)==Nparam*prop=1 : ",round(logLik(b_iubf_dr1),10))

## trigger situations where neighborhood is not contemporaneous
iubf(b_model5,Nubf=2, Nparam = 3,Nrep_per_param=3,nbhd=b_bag_nbhd_lookback1,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.8)
iubf(b_model5,Nubf=2, Nparam = 3,Nrep_per_param=3,nbhd=b_bag_nbhd_lookback2,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.8)

try(iubf(b_model5,Nubf=2, Nparam = 2,Nrep_per_param=3,nbhd=b_bag_nbhd,rw.sd=b_rw.sd,cooling.fraction.50=0.5,prop=0.8))

## --------------------------------------------
## using bm to test simulate and plot
## ____________________________________________

set.seed(9)

b_sim1 <- simulate(b_model,nsim=2,format='data.frame')
head(b_sim1,10)
b_sim2 <- simulate(b_model,nsim=2,format='data.frame',include.data=TRUE)
head(b_sim2,10)
b_sim3 <- simulate(b_model,nsim=2,format='spatPomps')

png(filename="bm-%02d.png",res=100)
plot(b_model,type="l",log=FALSE)
b_sim3v2 <- b_sim3[[1]]
b_sim3v2@data <- exp(b_sim3v2@data)
plot(b_sim3v2,type="l",log=TRUE)
plot(b_sim3[[2]],type="h",plot_unit_names=FALSE)
dev.off()
plot(b_sim3[[2]],type="h",plot_unit_names=TRUE) -> b_sim3_plot

print(b_model)

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
b_p_3d <- b_p
dim(b_p_3d) <- c(5,1,1)
dimnames(b_p_3d) <- c(dimnames(b_p),NULL)
vec_rmeasure(b_model,x=b_s,time=1, params=b_p_3d)

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
munit_measure(b_u, x=b_s, vc=b_vc, Np=1, unit = 2, time=1,
  params=b_array.params)
dunit_measure(b_u, y=b_y,x=b_s, unit=2, time=1, params=b_p)
runit_measure(b_u, x=b_s, unit=2, time=1, params=b_p)

## -------------------------------------------------------------
## test edge behaviors of basic components  
## _____________________________________________________________

dmeasure(b_model_zero_dmeasure,x=states(b_model))

dmeasure(b_model_zero_dmeasure,x=states(b_model),log=T)

vec_dmeasure(b_model_zero_dmeasure,y=obs(b_model_zero_dmeasure),
  x=states(b_model),units=1:U,
  times=1:2,params=coef(b_model_zero_dmeasure),log=T)[,,1]

## trigger error messages in dunit_measure.c
dunit_measure(b_model_no_dunit_measure, y=b_y,
  x=b_s, unit=1, time=1, params=b_p)
try(dunit_measure(b_model, y=b_y,
  x=b_s, unit=1, time=1:10, params=b_p))  
b_s2 <- 1:6
dim(b_s2) <- c(2,3,1)
b_s3 <- 1:6
dim(b_s3) <- c(1,3,2)
dimnames(b_s2) <- list(variable=dimnames(states(b_model))[[1]], rep=NULL)
b_p2 <- c(coef(b_model),coef(b_model))
dim(b_p2) <- c(length(coef(b_model)),2)
dimnames(b_p2) <- list(param=names(coef(b_model)))
try(dunit_measure(b_model, y=b_y, x=b_s2, unit=1, time=1, params=b_p2))
try(dunit_measure(b_model, y=b_y, x=b_s3, unit=1, time=1, params=b_p2))

## trigger error messages in runit_measure.c
runit_measure(b_model_no_runit_measure, x=b_s, unit=2, time=1, params=b_p)
try(runit_measure(b_model_no_runit_measure, x=b_s, unit=2, time=numeric(0), params=b_p))
try(runit_measure(b_model_no_runit_measure, x=b_s, unit=2, time=1:10, params=b_p))
try(runit_measure(b_model_no_runit_measure, x=b_s2, unit=2, time=1,params=b_p2))

## trigger error messages in vunit_measure.c
vunit_measure(b_model_no_vunit_measure, x=b_s, unit=2, time=1, params=b_p)
try(vunit_measure(b_model, x=b_s, unit=2, time=1:10, params=b_p))
try(vunit_measure(b_model, x=b_s2, unit=2, time=1, params=b_p2))
try(vunit_measure(b_model, x=b_s3, unit=2, time=1, params=b_p2))

## trigger error messages in eunit_measure.c
eunit_measure(b_model_no_eunit_measure, x=b_s, unit=2, time=1, params=b_p)
try(eunit_measure(b_model, x=b_s, unit=2, time=1:10, params=b_p))
try(eunit_measure(b_model, x=b_s2, unit=2, time=1, params=b_p2))

## trigger error messages in munit_measure.c
munit_measure(b_model_no_munit_measure, x=b_s, vc=b_vc, Np=1, unit = 2, time=1, params=b_array.params)
try(munit_measure(b_model, x=b_s, vc=b_vc, Np=1, unit = 2, time=1:10,params=b_array.params))
b_array.params2 <- array(c(b_p,b_p),
  dim = c(length(b_p),length(unit_names(b_model)), 2, 1),
  dimnames = list(params = rownames(b_p)))
try(munit_measure(b_model, x=b_s2, vc=b_vc, Np=3, unit = 2, time=1,params=b_array.params2))

## trigger error messages in fcstsampvar.c
try(.Call("do_fcst_samp_var",b_model,b_s,3,1:10,b_array.params,FALSE))
try(.Call("do_fcst_samp_var",b_model,b_s2,3,1,b_array.params2,FALSE))


## test spatPomp_Csnippet variable construction
spatPomp_Csnippet("lik=u;",unit_statenames="A",unit_obsnames=c("B","C"),
  unit_covarnames="D",
  unit_ivpnames="E",unit_paramnames="F",unit_vfnames="G")

## --------------------------------------------
## using bm to test spatPomp() replacement functionality
## ____________________________________________

b_rep1 <- spatPomp(b_model,params=coef(b_model))
for(slt in slotNames(b_model)) if(!identical(slot(b_model,slt),
  slot(b_rep1,slt))) print(slt)

# test parameter replacement
b_rep2 <- spatPomp(b_model,params=coef(b_model)+1)
if(!identical(coef(b_rep2),coef(b_model)+1)) stop('problem with parameter replacement')

# test do-nothing behavior
b_rep3 <- spatPomp(b_model)

## --------------------------------------------
## using bm to test spatPomp() warning messages
## ____________________________________________

print("The following deliver error messages, to test them")
try(spatPomp(data=as.data.frame(b_model),units=NULL),outFile=stdout())
try(spatPomp("test on type character"))

try(spatPomp())
b_data <- as.data.frame(b_model)
try(spatPomp(data=b_data,times="time",units="unit"))
try(spatPomp(data=b_data,times="NONSENSE",units="unit",t0=0))
try(spatPomp(data=b_data,times="time",units="NONSENSE",t0=0))
spatPomp(data=b_data,times="time",units="unit",t0=0,
  params=list(coef(b_model)))
b_data2 <- b_data
names(b_data2) <- c("time","unit","X","X")
try(spatPomp(data=b_data2,times="time",units="unit"))
b_data_only_model <- spatPomp(data=b_data,times="time",units="unit",
  t0=0)

# test error messages for covariates with data.frame class for spatPomp()
b_covar_error <- data.frame(time_name_error=0:2,Z=3:5)
try(spatPomp(data=b_data,times="time",units="unit",t0=0,
  covar=b_covar_error))

b_unit_covar_names_error <- data.frame(time=c(0:2,0:2),
  JUNK=rep(c("U1","U2"),each=3), Z=rep(3:5,times=2))
try(spatPomp(data=b_data,times="time",units="unit",t0=0,
  covar=b_unit_covar_names_error))

b_shared_covar <- data.frame(time=0:2,Z=3:5)
model_shared_covar <- spatPomp(data=b_data,times="time",units="unit",
  t0=0,covar=b_shared_covar, shared_covarnames="Z")
try(as.data.frame(model_shared_covar))

b_unit_covar <- data.frame(time=c(0:2,0:2),unit=rep(c("U1","U2"),each=3),
  Z=rep(3:5,times=2))
model_unit_covar <- spatPomp(data=b_data,times="time",units="unit",
  t0=0,covar=b_unit_covar,skeleton=NULL,partrans=NULL,
  unit_accumvars = "JUNK")

try(spatPomp(data=b_data,times="time",units="unit",t0=0,
  covar=b_unit_covar,shared_covarnames ="JUNK"))

# test spatPomp warnings with argument of class spatPomp

# perhaps surprisingly, this gives no error
spatPomp(model_unit_covar,timename="JUNK",unitname="JUNK",
  unit_accumvars="JUNK", globals=Csnippet("JUNK"),
  partrans=NULL,skeleton=NULL)
  
try(spatPomp(data=model_unit_covar,covar=b_covar_error))

spatPomp(model_shared_covar)

spatPomp(data=model_shared_covar,covar=b_shared_covar,
  shared_covarnames="Z")

spatPomp(data=model_unit_covar,covar=b_unit_covar)

try(spatPomp(data=model_unit_covar,covar=b_shared_covar))

try(spatPomp(data=model_unit_covar,covar=b_unit_covar,
  shared_covarnames="Z"))


## --------------------------------------------
## test utility functions
## ____________________________________________

spatPomp:::undefined()
spatPomp:::undefined(NULL)
spatPomp:::undefined("JUNK")

.Call("spatPomp_systematic_resampling",c(0.1,0.2),2)
try(.Call("spatPomp_systematic_resampling",c(-0.1,-0.2),2))

try(spatPomp:::conc())
try(spatPomp:::conc("a","b"))

(t_spatPompList <- is(c(b_model, b_model), "spatPompList"))

class( c(c(b_model,b_model),b_model))

(t_bpfilterList <- is(
  c(b_bpfilter,b_bpfilter),
  "bpfilterList"
))

## ibpfilterList is tested in he10

stopifnot(all(t_spatPompList, t_bpfilterList))



