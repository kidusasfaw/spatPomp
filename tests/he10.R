library(spatPomp)
set.seed(22)

#model_type <- "he10"
model_type <- "mostly shared"
parNames <- c("alpha","R0","g","sigma","gamma","amplitude","cohort","sigmaSE","S_0","E_0","I_0","rho","psi","iota","mu")
# he10 defaults to alpha=1, cohort=0, which means the usual transformations are undefined.
# here, we don't estimate either
if(model_type == "mostly fixed"){
  sharedParNames <- c("R0","psi")
  unitParNames <- c("rho","S_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
} else if(model_type == "mostly shared"){
  sharedParNames <- c("R0","psi","g","sigma","gamma","amplitude","sigmaSE")
  unitParNames <- c("rho","S_0","E_0","I_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
} else if(model_type == "plausible parameters shared"){
  # parameters are shared when that makes mechanistic sense.
  sharedParNames <- c("R0","g","sigma","gamma","amplitude")
  unitParNames <- c("sigmaSE","S_0","E_0","I_0","rho","psi")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
} else if(model_type == "all unit-specific"){
  # all parameters estimated except life expecancy
  # and immigration, which should not be needed when there is coupling
  fixedParNames <- c("mu","iota")
  sharedParNames <- NULL
  unitParNames <- setdiff(parNames,fixedParNames)
  estParNames <- c(sharedParNames,unitParNames)
} else if(model_type == "he10"){
  # all the parameters estimated by He et al (2010) Table 2
  fixedParNames <- c("mu","g")
  sharedParNames <- NULL
  unitParNames <- setdiff(parNames,fixedParNames)
  estParNames <- c(sharedParNames,unitParNames)
}

## test error messages
try(he10(U=5,towns_selected="JUNK"))
try(he10(U=1000,towns_selected=1:1000))
try(he10(U=5,Tmax=2024))

## Note: here we assume that there are no unestimated unit-specific
## parameters. That could readily be accommodated if needed.

h_model <- he10(U=2,dt=4/365,Tmax=1950.1,
  expandedParNames=estParNames)

coef(h_model)

h_bpfilter <- bpfilter(h_model,Np=10,block_size=1)

paste("bpfilter logLik for he10 model:",logLik(h_bpfilter))


h_U <- length(unit_names(h_model))

ivpParNames <- c("S_0","E_0","I_0")
ivpEstParNames <- intersect(ivpParNames,estParNames)
regEstParNames <- setdiff(estParNames,ivpParNames)

estParNames_expanded <- unlist(lapply(estParNames,function(x)paste0(x,1:h_U)))
regEstParNames_expanded <- unlist(lapply(regEstParNames,function(x)paste0(x,1:h_U)))
ivpEstParNames_expanded <- unlist(lapply(ivpEstParNames,function(x)paste0(x,1:h_U)))
fixedParNames_expanded <- paste0(fixedParNames,1)


reg_rw.sd <- rep(list(0.02),times=length(regEstParNames_expanded))
names(reg_rw.sd) <- regEstParNames_expanded
if("alpha"%in%estParNames) reg_rw.sd[paste0("alpha",1:h_U)] <- 0.005

ivp_rw.sd <- lapply(ivpEstParNames_expanded,function(x)expression(ivp(0.05)))
names(ivp_rw.sd) <- ivpEstParNames_expanded
h_rw.sd <- do.call(rw_sd,c(reg_rw.sd,ivp_rw.sd))

all_units = seq_len(length(unit_names(h_model)))
nblocks = 2
block_list = split(all_units, sort(all_units %% nblocks))
block_list <- lapply(block_list, as.integer)

set.seed(3)
h_ibpf <- ibpf(h_model,
  params=coef(h_model),
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  Nbpf=2,
  spat_regression=0.1,
  Np=10,
  rw.sd=h_rw.sd,
  cooling.fraction.50=0.5,
  block_list=block_list
)

h_bpfilter <- bpfilter(h_ibpf,Np=10,block_size=1)

paste("ibpf logLik for he10 model:",logLik(h_bpfilter))

# test whether specifying Np as a function gives the same result
set.seed(3)
h_ibpf2 <- ibpf(
h_model,
  params=coef(h_model),
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  Nbpf=2,
  spat_regression=0.1,
  Np=function(k) 10,
  rw.sd=h_rw.sd,
  cooling.fraction.50=0.5,
  block_list=block_list
)

h_bpfilter2 <- bpfilter(h_ibpf2,Np=10,block_size=1)

if (logLik(h_bpfilter2)!=logLik(h_bpfilter))
  stop("in ibpf: Np specified as a function gives a different result from Np as a scalar")
  
coef(h_ibpf)

# test errors for ibpf on class 'missing' or character
try(ibpf())
try(ibpf("h_model"))

# test errors for ibpf on class spatPomp
try(ibpf(h_model))
try(ibpf(h_model,Nbpf=2))
try(ibpf(h_model,Nbpf=2,rw.sd=rw_sd(mu1=0.1)))
try(ibpf(h_model,Nbpf=NA,Np=10))
try(ibpf(h_model,Nbpf=NA,Np=10,block_size=1))
try(ibpf(h_model,Nbpf=NA,Np=10,block_size=1,sharedParNames=NULL))
try(ibpf(h_model,Nbpf=2,rw.sd=rw_sd(mu1=0.1),Np=10,sharedParNames=sharedParNames,
  unitParNames=unitParNames))
try(ibpf(h_model,Nbpf=2,rw.sd=rw_sd(mu1=0.1),Np=10,sharedParNames=sharedParNames,
  unitParNames=unitParNames,block_list=block_list,block_size=1))
try(ibpf(h_model,Nbpf=2,rw.sd=rw_sd(mu1=0.1),Np=10,sharedParNames=sharedParNames,
  unitParNames=unitParNames,block_list=block_list))
try(ibpf(h_model,Nbpf=2,rw.sd=rw_sd(mu1=0.1),Np=5,sharedParNames=sharedParNames,
  unitParNames=unitParNames,spat_regression=0.5,block_size=10))

try(ibpf(h_model,Nbpf=NULL,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.1)))
try(ibpf(h_model,Nbpf=NULL,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.1),
  sharedParNames=NULL))
try(ibpf(h_model,Nbpf=NULL,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.1),
  sharedParNames=NULL,unitParNames=NULL))
try(ibpf(h_model,Nbpf=NULL,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.1),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=0.5))
try(ibpf(h_model,Nbpf=NULL,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.00001),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=0.5,
  spat_regression=0.5))
try(ibpf(h_model,Nbpf=1,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.00001),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=12,
  spat_regression=0.5))
try(ibpf(h_model,Nbpf=-1,block_list=block_list,Np=10,rw.sd=rw_sd(mu1=0.00001),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=0.5,
  spat_regression=0.5))

# test errors on Np specification
try(ibpf(h_model,Nbpf=2,block_list=block_list,Np=NULL,rw.sd=rw_sd(mu1=0.00001),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=0.5,
  spat_regression=0.5))
try(ibpf(h_model,Nbpf=2,block_list=block_list,Np=1:100,rw.sd=rw_sd(mu1=0.00001),
  sharedParNames=sharedParNames,unitParNames=unitParNames,cooling.fraction.50=0.5,
  spat_regression=0.5))
try(ibpf(h_model,Nbpf=2,block_list=block_list,Np="a character vector",
  rw.sd=rw_sd(mu1=0.00001),sharedParNames=sharedParNames,
  unitParNames=unitParNames,cooling.fraction.50=0.5,spat_regression=0.5))
try(ibpf(h_model,Nbpf=2,block_list=block_list,Np=c(10,10),
  rw.sd=rw_sd(mu1=0.00001),sharedParNames=sharedParNames,
  unitParNames=unitParNames,cooling.fraction.50=0.5,spat_regression=0.5))

# test ibpf errors on class ibpfd_spatPomp

capture.output(ibpf(h_ibpf,sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  .paramMatrix=h_ibpf@paramMatrix,verbose=TRUE)) -> out
try(ibpf(h_ibpf,block_size="JUNK",block_list="JUNK"))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  block_size=1,Nbpf <- 0.1))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  block_size=3))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  Np=function(n) "JUNK"))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  Np=function(n) -1))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  .paramMatrix=h_ibpf@paramMatrix,Np=7))
try(ibpf(h_ibpf,sharedParNames=sharedParNames,unitParNames=unitParNames,
  .paramMatrix=h_ibpf@paramMatrix[,1,drop=FALSE],Np=1))

# test ibpf on class bpfilterd_spatPomp
try(ibpf(h_bpfilter,block_list=block_list,block_size=1))
try(ibpf(h_bpfilter,block_size=23))
try(ibpf(h_bpfilter))


# test ibpf with missing basic model component
h_model2 <- spatPomp(h_model,rprocess=NULL)
try(h_ibpf2 <- ibpf(h_model2,
  params=coef(h_model),
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  Nbpf=2,
  spat_regression=0.1,
  Np=10,
  rw.sd=h_rw.sd,
  cooling.fraction.50=0.5,
  block_list=block_list
))

## test error message when munit_measure is undefined
try(girf(h_model,kind="moment",
  Np=10,Ninter=2,Nguide=10,lookahead=1,tol=1e-5))

# Create second ibpfd_spatPomp object with different chain length, to test error
h_ibpf3 <- ibpf(h_model,
                params=coef(h_model),
                sharedParNames=sharedParNames,
                unitParNames=unitParNames,
                Nbpf=3,
                spat_regression=0.1,
                Np=10,
                rw.sd=h_rw.sd,
                cooling.fraction.50=0.5,
                block_list=block_list
)

# Should correctly make ibpfList object
is(c(h_ibpf, h_ibpf), "ibpfList")

# Throws error because they have different chain lengths
try(c(h_ibpf, h_ibpf3))

