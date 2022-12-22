library(spatPomp)
set.seed(0)

b2_U <- 4
b2 <- bm2(U=b2_U,N=2,unit_specific_names="rho")

paste("bpfilter logLik for bm2 model:",logLik(bpfilter(b2,Np=5,block_size=1)))

# here there are no transformations so use small rw.sd. to avoid negatives
b2_rw_list <- rep(list(0.001),times=b2_U) 
names(b2_rw_list) <-paste0("rho",1:b2_U)
b2_rw.sd <- do.call(rw.sd,b2_rw_list)

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
paste("check ibpf on ipfd_spatPomp: ",
  logLik(b2_ibpf)==logLik(b2_ibpf_repeat))



