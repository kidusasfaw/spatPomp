library(spatPomp)

## test error message
try(bm2(U=3,N=2,unit_specific_names="rho",shared_names="sigma"))
try(bm2_kalman_logLik(bm2(U=1,N=2)))

b2_U <- 4

set.seed(0)
b2 <- bm2(U=b2_U,N=2,unit_specific_names="rho")

set.seed(0)
b2alt <- bm2(U=b2_U,N=2,shared_names=c("tau","sigma","X_0"))
if(any(obs(b2)!=obs(b2alt))) stop("two different ways to specify the same bm2 model should match")

b2_bpfilter <- bpfilter(b2,Np=5,block_size=1)
paste("bpfilter logLik for bm2 model:",logLik(b2_bpfilter))

# here there are no transformations so use small rw.sd. to avoid negatives
b2_rw_list <- rep(list(0.001),times=b2_U) 
names(b2_rw_list) <-paste0("rho",1:b2_U)
b2_rw.sd <- do.call(rw_sd,b2_rw_list)

b2_units = seq_len(b2_U)
b2_nblocks = b2_U/2
b2_block_list = split(b2_units, sort(b2_units %% b2_nblocks))
b2_block_list <- lapply(b2_block_list, as.integer)

b2_ibpf <- ibpf(b2,
  params=coef(b2),
  sharedParNames=NULL,
  unitParNames="rho",
  Nbpf=2,
  spat_regression=0.1,
  Np=5,
  rw.sd=b2_rw.sd,
  cooling.fraction.50=0.5,
  block_list=b2_block_list
)

paste("ibpf logLik for b2 model:",logLik(b2_ibpf))

paste("kf logLik for b2:",bm2_kalman_logLik(b2))

#######################################################################
# test ibpf with argument re-use and replacement
# check that results match after re-use and replacement
#######################################################################

set.seed(5)
b2_ibpf <- ibpf(b2,Np=5,block_size=1,Nbpf=2,
  rw.sd=b2_rw.sd,
  cooling.frac=0.5, spat_regression=0.1,
  unitParNames="rho",sharedParNames=NULL)

paste("bm2 ibpf loglik: ",round(logLik(b2_ibpf),10))

set.seed(5)
b2_ibpf_repeat <- ibpf(b2_ibpf,params=coef(b2), unitParNames="rho",
  sharedParNames=NULL,spat_regression=0.1)
paste("check ibpf on ibpfd_spatPomp: ",
  logLik(b2_ibpf_repeat)==logLik(b2_ibpf))

set.seed(5)
b2_ibpf_bpfilterd <- ibpf(b2_bpfilter,Np=5,block_size=1,Nbpf=2,
  rw.sd=b2_rw.sd,
  cooling.frac=0.5, spat_regression=0.1,
  unitParNames="rho",sharedParNames=NULL)
paste("check ibpf on bpfilterd_spatPomp: ",
  logLik(b2_ibpf_bpfilterd)==logLik(b2_ibpf))



