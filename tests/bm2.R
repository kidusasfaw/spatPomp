library(spatPomp)
set.seed(0)

b2_U <- 4
b2 <- bm2(U=b2_U,N=5,unit_specific_names="rho")

paste("bpfilter logLik for bm2 model:",logLik(bpfilter(b2,Np=10,block_size=1)))

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
  Np=10,
  rw.sd=b2_rw.sd,
  cooling.fraction.50=0.5,
  block_list=b2_block_list
)

paste("ibpf logLik for b2 model:",logLik(bpfilter(b2_ibpf,Np=10,block_size=1)))





