library(spatPomp)
set.seed(22)

model_type <- "he10"
parNames <- c("alpha","R0","g","sigma","gamma","amplitude","cohort","sigmaSE","S_0","E_0","I_0","rho","psi","iota","mu")
if(model_type == "mostly fixed"){
  sharedParNames <- c("R0","psi")
  unitParNames <- c("rho","S_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
} else if(model_type == "mostly shared"){
  sharedParNames <- c("alpha","R0","psi","g","sigma","gamma","amplitude","cohort","sigmaSE")
  unitParNames <- c("rho","S_0","E_0","I_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
} else if(model_type == "plausible parameters shared"){
  # parameters are shared when that makes mechanistic sense. 
  sharedParNames <- c("alpha","R0","g","sigma","gamma","amplitude","cohort")
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

# Note: here we assume that there are no unestimated unit-specific
# parameters. That could readily be accommodated if needed.

h_model <- he10(,U=2,dt=4/365,Tmax=1950.5,
  expandedParNames=estParNames)

coef(h_model)

paste("bpfilter logLik for he10 model:",logLik(bpfilter(h_model,Np=10,block_size=1)))


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
h_rw.sd <- do.call(rw.sd,c(reg_rw.sd,ivp_rw.sd))

all_units = seq_len(length(unit_names(h_model)))
nblocks = 2
block_list = split(all_units, sort(all_units %% nblocks))
block_list <- lapply(block_list, as.integer)

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

paste("ibpfilter logLik for he10 model:",logLik(bpfilter(h_model,Np=10,block_size=1)))





