## ----cores,echo=F,cache=F,message = FALSE, warnings = FALSE,results='hide'----
library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
ggplot2::theme_set(ggplot2::theme_bw())


## ----packages,include=F,echo=F,cache=F----------------------------------------
library("spatPomp")
library("ggplot2")
library("tidyverse")
library("knitr")
stopifnot(packageVersion("pomp")>="5.0")






## ----setup,cache=F,include=F--------------------------------------------------
myround<- function (x, digits = 1) {
  # adapted from the broman package
  # solves the bug that round() kills significant trailing zeros
  if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  if (digits < 1) {
    as.character(round(x,digits))
  } else {
    tmp <- sprintf(paste("%.", digits, "f", sep = ""), x)
    zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
    tmp[tmp == paste0("-", zero)] <- zero
    tmp
  }
}
mysignif <- function (x, digits = 1) {
  myround(x, digits - ceiling(log10(abs(x))))
}


## ----bm-setup,echo=F,eval=T---------------------------------------------------
set.seed(20)


## ----bm-model,echo=T,eval=T---------------------------------------------------
library(spatPomp)
i <- 2
U <- 4
b <- bm2(U=U,N=switch(i,10,200),unit_specific_names="rho")
plot(b)


## ----bm-followup,echo=F,eval=T------------------------------------------------
if(i==1) cores <- 2   # no need to use all available cores when i=1




## ----bm_logLik_eval,cache=F,echo=F,eval=T-------------------------------------
bm_dir <- paste0("bm_",i,"/")
if(!dir.exists(bm_dir)) dir.create(bm_dir)
stew(file=paste0(bm_dir,"logLik.rda"),seed=10,{

  cat(capture.output(sessionInfo()),
    file=paste0(bm_dir,"sessionInfo.txt"),sep="\n")

kf_logLik <- bm2_kalman_logLik(b)
pf_logLik <- replicate(10,
  logLik(pfilter(b,switch(i,10,1000)))
)
bpf_logLik2 <- replicate(10,
  logLik(bpfilter(b,switch(i,10,1000),block_size=2))
)
  bpf_logLik1 <- replicate(10,
    logLik(bpfilter(b,switch(i,10,1000),block_size=1))
  )
  bpf_logLik4 <- replicate(10,
    logLik(bpfilter(b,switch(i,10,1000),block_size=4))
  )
})


## ----kf_mle,echo=F,eval=T-----------------------------------------------------
stew(file=paste0(bm_dir,"kf_mle.rda"),seed=999,{
  bm2_negLogLik <- function(rho){
    coef(b,names(rho)) <- unname(rho)
    -bm2_kalman_logLik(b)
  }

  rho_init <- c(rho1=0.4,rho2=0.4,rho3=0.4,rho4=0.4)
  bm2_negLogLik(rho_init)
  bm2_mle <- optim(rho_init,bm2_negLogLik)
})


## ----kf_prof,echo=F,eval=T----------------------------------------------------
stew(file=paste0(bm_dir,"kf_prof.rda"),seed=999,{
bm2_negProf <- function(rho_est,rho_fixed,u){
  rho_new <- rep(NA,U)
  names(rho_new) <- paste0("rho",1:U)
  rho_new[names(rho_est)] <- unname(rho_est)
  rho_new[u] <- rho_fixed
  coef(b,names(rho_new)) <- unname(rho_new)
  -bm2_kalman_logLik(b)
}

prof_init <- rep(0.4,U)
names(prof_init) <- paste0("rho",1:U)
prof_vals <- seq(from=0,to=0.8,length.out=switch(i,5,20))

kf_prof <- function(rho_vals,u){
  foreach(rho_prof=rho_vals)%dopar% {
    optim(prof_init[-u],bm2_negProf,rho_fixed=rho_prof,u=u)
  }
}

prof1 <- kf_prof(prof_vals,u=1)
prof2 <- kf_prof(prof_vals,u=2)
prof3 <- kf_prof(prof_vals,u=3)
prof4 <- kf_prof(prof_vals,u=4)

prof1_logLik <- sapply(prof1,function(x)-x$value)
prof2_logLik <- sapply(prof2,function(x)-x$value)
prof3_logLik <- sapply(prof3,function(x)-x$value)
prof4_logLik <- sapply(prof4,function(x)-x$value)
})

prof_range <- range(c(-bm2_mle$value,-bm2_mle$value-10))
if(0){
  par(mfrow=c(2,2))
  plot(y=prof1_logLik,x=prof_vals,ty="l",ylim=prof_range,
    ylab="Log-likelihood",xlab="rho1")
  plot(y=prof2_logLik,x=prof_vals,ty="l",ylim=prof_range,
    ylab="Log-likelihood",xlab="rho2")
  plot(y=prof3_logLik,x=prof_vals,ty="l",ylim=prof_range,
    ylab="Log-likelihood",xlab="rho3")
  plot(y=prof4_logLik,x=prof_vals,ty="l",ylim=prof_range,
    ylab="Log-likelihood",xlab="rho4")
}




## ----ibpf-mle-eval,eval=T,echo=F----------------------------------------------
stew(file=paste0(bm_dir,"ibpf_mle.rda"),seed=999,{
rho_start <- seq(from=0.2,to=0.8,length.out=U)
params_start <- coef(b)
params_start[paste0("rho",1:U)] <- rho_start
ibpf_mle_searches <- foreach(reps=1:switch(i,3,10))%dopar%{
  ibpf(b,params=params_start,
    Nbpf=switch(i,2,50),Np=switch(i,10,1000),
    rw.sd=rw_sd(rho1=0.02,rho2=0.02,rho3=0.02,rho4=0.02),
    unitParNames="rho",
    sharedParNames=NULL,
    block_size=2,
    cooling.fraction.50=0.5
  )
}
})


## ----ibpf-mle-lik,eval=T,echo=F-----------------------------------------------
stew(file=paste0(bm_dir,"ibpf_mle_lik.rda"),seed=878,{
ibpf_kf_eval <- sapply(ibpf_mle_searches,bm2_kalman_logLik)
ibpf_bpf_eval <- foreach(bm2fit=ibpf_mle_searches,.combine=cbind) %dopar% {
  replicate(switch(i,3,10),
    logLik(bpfilter(bm2fit,block_size=2,Np=switch(i,10,1000)))
  )
}
})


## ----ibpf-prof,eval=T,echo=F--------------------------------------------------
stew(file=paste0(bm_dir,"prof1.rda"),seed=722,{
  tic <- Sys.time()
  u <- 1
  rho_start <- seq(from=0.2,to=0.8,length.out=U)
  reps <- switch(i,3,10)
  param_start_matrix <- matrix(params_start,byrow=TRUE,
    ncol=length(coef(b)),nrow=reps*length(prof_vals),
    dimnames=list(NULL,names(coef(b))))
  param_start_matrix[,paste0("rho",u)] <- rep(prof_vals,each=reps)
  ibpf_prof_searches <- foreach(s=1:nrow(param_start_matrix))%dopar%{
    rw_sd_call <- rep(list(0.02),times=U-1)
    names(rw_sd_call) <- lapply((1:U)[-u],function(x)paste0("rho",x))
    ibpf(b,
      params=param_start_matrix[s,],
      Nbpf=switch(i,2,50),Np=switch(i,20,1000),
      rw.sd=do.call(rw_sd,rw_sd_call),
      unitParNames="rho",
      sharedParNames=NULL,
      block_size=2,
      cooling.fraction.50=0.5
    )
  }
  toc <- Sys.time()
})
prof1time <- toc-tic


## ----ibpf-prof-lik,eval=T,echo=F----------------------------------------------
stew(file=paste0(bm_dir,"prof1_lik.rda"),seed=868,{
tic <- Sys.time()
prof1_kf_eval <- sapply(ibpf_prof_searches,bm2_kalman_logLik)
ibpf_bpf_eval <- foreach(bm2fit=ibpf_prof_searches,.combine=cbind) %dopar% {
  replicate(switch(i,5,20),
    logLik(bpfilter(bm2fit,block_size=2,Np=switch(i,10,1000)))
  )
}
ibpf_pf_eval <- foreach(bm2fit=ibpf_prof_searches,.combine=cbind) %dopar% {
  replicate(switch(i,5,20),
    logLik(pfilter(bm2fit,Np=switch(i,10,1000)))
  )
}
toc <- Sys.time()
prof1cores <- cores
})
prof1eval_time <- toc-tic


## ----ibpf-prof-rho-plot, fig.height=8, fig.width=6, out.width="5.5in", fig.cap = "Top: profile for $\\rho_1$ using an IBPF search with likelihood computed using BPF, for $K=2$ blocks each having 2 units. Middle: Exact profile (dashed red line) and the same IBPF search with likelihood computed exactly using the Kalman filter. Bottom: The same IBPF search with likelihood computed using the particle filter. Vertical lines show the MLE and a 95\\% confidence interval, with a dotted line at the true parameter value.",eval=T,echo=F----

par(mfrow=c(3,1))
par(mai=c(0.9,0.9,0.1,0.2))

##### bpf evaluation
p1bpf <- apply(matrix(apply(ibpf_bpf_eval,2,mean),nrow=reps),2,max)
bpf_prof_range <- c(max(p1bpf)-10,max(p1bpf))
plot(y=p1bpf,x=prof_vals,ylim=bpf_prof_range,
  ylab="Log-likelihood",xlab="")
p1bpf_mcap <- mcap(logLik=p1bpf,parameter=prof_vals)
lines(p1bpf_mcap$fit$parameter,p1bpf_mcap$fit$smoothed)

abline(v=p1bpf_mcap$ci)
abline(v=p1bpf_mcap$mle)
abline(v=coef(b)["rho1"],lty="dotted")

###### kf evaluation
p1kf <- apply(matrix(prof1_kf_eval,nrow=reps),2,max)
plot(y=p1kf,x=prof_vals,ylim=prof_range,
  ylab="Log-likelihood",xlab="")
p1kf_mcap <- mcap(logLik=p1kf,parameter=prof_vals)
lines(p1kf_mcap$fit$parameter,p1kf_mcap$fit$smoothed)
p1exact_mcap <- mcap(logLik=prof1_logLik,prof_vals)
# points(y=prof1_logLik,x=prof_vals,col="red")
lines(p1exact_mcap$fit$parameter,p1exact_mcap$fit$smoothed,col="red",lty="dashed")

abline(v=p1kf_mcap$ci)
abline(v=p1exact_mcap$ci,col="red",lty="dashed")
abline(v=p1kf_mcap$mle)
abline(v=coef(b)["rho1"],lty="dotted")

######## pf evaluation
p1pf <- apply(matrix(apply(ibpf_pf_eval,2,mean),nrow=reps),2,max)
pf_prof_range <- c(max(p1pf)-10,max(p1pf))
plot(y=p1pf,x=prof_vals,ylim=pf_prof_range,
  ylab="Log-likelihood",xlab="rho1")
p1pf_mcap <- mcap(logLik=p1pf,parameter=prof_vals)
lines(p1pf_mcap$fit$parameter,p1pf_mcap$fit$smoothed)

abline(v=p1pf_mcap$ci)
abline(v=p1pf_mcap$mle)
abline(v=coef(b)["rho1"],lty="dotted")



## ----bm-sigma-model,echo=F,eval=T---------------------------------------------
set.seed(20)
b_sig <- bm2(U=4,N=switch(i,10,200),unit_specific_names="sigma")

bm2_sig_negLogLik <- function(sigma){
    coef(b_sig,names(sigma)) <- unname(sigma)
    -bm2_kalman_logLik(b_sig)
  }

stew(file=paste0(bm_dir,"kf_mle_sig.rda"),seed=256,{
  sigma_init <- c(sigma1=1,sigma2=1,sigma3=1,sigma4=1)
  bm2_sig_negLogLik(sigma_init)
  bm2_sig_mle <- optim(sigma_init,bm2_sig_negLogLik)
})



## ----kf-prof-sig,echo=F,eval=T------------------------------------------------
stew(file=paste0(bm_dir,"kf_prof_sig.rda"),seed=512,{

bm2_sig_negProf <- function(sigma_est,sigma_fixed,u){
  sigma_new <- rep(NA,U)
  names(sigma_new) <- paste0("sigma",1:U)
  sigma_new[names(sigma_est)] <- unname(sigma_est)
  sigma_new[u] <- sigma_fixed
  coef(b_sig,names(sigma_new)) <- unname(sigma_new)
  -bm2_kalman_logLik(b_sig)
}

prof_sig_init <- rep(1,U)
names(prof_sig_init) <- paste0("sigma",1:U)
prof_sig_vals <- seq(from=0.4,to=1.6,length.out=switch(i,5,20))

kf_sig_prof <- function(sigma_vals,u){
  foreach(sigma_prof=sigma_vals)%dopar% {
    optim(prof_sig_init[-u],bm2_sig_negProf,sigma_fixed=sigma_prof,u=u)
  }
}

prof_sig1 <- kf_sig_prof(prof_sig_vals,u=1)
prof_sig2 <- kf_sig_prof(prof_sig_vals,u=2)
prof_sig3 <- kf_sig_prof(prof_sig_vals,u=3)
prof_sig4 <- kf_sig_prof(prof_sig_vals,u=4)

prof_sig1_logLik <- sapply(prof_sig1,function(x)-x$value)
prof_sig2_logLik <- sapply(prof_sig2,function(x)-x$value)
prof_sig3_logLik <- sapply(prof_sig3,function(x)-x$value)
prof_sig4_logLik <- sapply(prof_sig4,function(x)-x$value)
})

prof_sig_range <- range(c(-bm2_sig_mle$value,-bm2_sig_mle$value-10))
if(0){
  par(mfrow=c(2,2))
  plot(y=prof_sig1_logLik,x=prof_sig_vals,ty="l",ylim=prof_sig_range,
    ylab="Log-likelihood",xlab="sigma1")
  plot(y=prof_sig2_logLik,x=prof_sig_vals,ty="l",ylim=prof_sig_range,
    ylab="Log-likelihood",xlab="sigma2")
  plot(y=prof_sig3_logLik,x=prof_sig_vals,ty="l",ylim=prof_sig_range,
    ylab="Log-likelihood",xlab="sigma3")
  plot(y=prof_sig4_logLik,x=prof_sig_vals,ty="l",ylim=prof_sig_range,
    ylab="Log-likelihood",xlab="sigma4")
}


## ----ibpf-prof-sig,eval=T,echo=F----------------------------------------------
stew(file=paste0(bm_dir,"prof_sig1.rda"),seed=722,{
tic <- Sys.time()
u <- 1
sigma_start <- seq(from=0.4,to=1.6,length.out=U)
params_sig_start <- coef(b_sig)
params_sig_start[paste0("sigma",1:U)] <- sigma_start
reps <- switch(i,2,10)
param_sig_start_matrix <- matrix(params_sig_start,byrow=TRUE,
  ncol=length(coef(b_sig)),nrow=reps*length(prof_sig_vals),
  dimnames=list(NULL,names(coef(b_sig))))
param_sig_start_matrix[,paste0("sigma",u)] <- rep(prof_sig_vals,each=reps)
ibpf_prof_sig_searches <- foreach(s=1:nrow(param_sig_start_matrix))%dopar%{
  rw_sd_call <- rep(list(0.02),times=U-1)
  names(rw_sd_call) <- lapply((1:U)[-u],function(x)paste0("sigma",x))
  ibpf(b_sig,
    params=param_sig_start_matrix[s,],
    Nbpf=switch(i,2,50),Np=switch(i,10,1000),
    rw.sd=do.call(rw_sd,rw_sd_call),
    unitParNames="sigma",
    sharedParNames=NULL,
    block_size=2,
    cooling.fraction.50=0.5
  )
}
toc <- Sys.time()
})
prof_sig1time <- toc-tic


## ----ibpf-prof-sig-lik,eval=T,echo=F------------------------------------------
stew(file=paste0(bm_dir,"prof_sig1_lik.rda"),seed=868,{
tic <- Sys.time()
prof_sig1_kf_eval <- sapply(ibpf_prof_sig_searches,bm2_kalman_logLik)
ibpf_bpf_sig_eval <- foreach(bm2fit=ibpf_prof_sig_searches,.combine=cbind) %dopar% {
  replicate(switch(i,2,10),
    logLik(bpfilter(bm2fit,block_size=2,Np=switch(i,10,1000)))
  )
}
ibpf_pf_sig_eval <- foreach(bm2fit=ibpf_prof_sig_searches,.combine=cbind) %dopar% {
  replicate(switch(i,2,10),
    logLik(pfilter(bm2fit,Np=switch(i,10,1000)))
  )
}
toc <- Sys.time()
prof_sig_cores <- cores
})
prof_sig1eval_time <- toc-tic


## ----ibpf-prof-sig-plot, fig.height=8, fig.width=6, out.width="5.5in", fig.cap = "Top: profile for $\\sigma_1$ using an IBPF search with likelihood computed using BPF, for $K=2$ blocks each having 2 units. Middle: Exact profile (dashed red line) and the same IBPF search with likelihood computed exactly using the Kalman filter. Bottom: The same IBPF search with likelihood computed using the particle filter. Vertical lines show the MLE and a 95\\% confidence interval, with a dotted line at the true parameter value.",eval=T,echo=F----

par(mfrow=c(3,1))
par(mai=c(0.9,0.9,0.1,0.2))

##### bpf evaluation
p1_sig_bpf <- apply(matrix(apply(ibpf_bpf_sig_eval,2,mean),nrow=reps),2,max)
bpf_prof_sig_range <- c(max(p1_sig_bpf)-10,max(p1_sig_bpf))
plot(y=p1_sig_bpf,x=prof_sig_vals,ylim=bpf_prof_sig_range,
  ylab="Log-likelihood",xlab="")
p1_sig_bpf_mcap <- mcap(logLik=p1_sig_bpf,parameter=prof_sig_vals)
lines(p1_sig_bpf_mcap$fit$parameter,p1_sig_bpf_mcap$fit$smoothed)

abline(v=p1_sig_bpf_mcap$ci)
abline(v=p1_sig_bpf_mcap$mle)
abline(v=coef(b_sig)["sigma1"],lty="dotted")

###### kf evaluation
p1_sig_kf <- apply(matrix(prof_sig1_kf_eval,nrow=reps),2,max)
plot(y=p1_sig_kf,x=prof_sig_vals,ylim=prof_sig_range,
  ylab="Log-likelihood",xlab="")
p1_sig_kf_mcap <- mcap(logLik=p1_sig_kf,parameter=prof_sig_vals)
lines(p1_sig_kf_mcap$fit$parameter,p1_sig_kf_mcap$fit$smoothed)
p1_sig_exact_mcap <- mcap(logLik=prof_sig1_logLik,prof_sig_vals)
# points(y=prof_sig1_logLik,x=prof_sig_vals,col="red")
lines(p1_sig_exact_mcap$fit$parameter,p1_sig_exact_mcap$fit$smoothed,col="red",lty="dashed")

abline(v=p1_sig_kf_mcap$ci)
abline(v=p1_sig_exact_mcap$ci,col="red",lty="dashed")
abline(v=p1_sig_kf_mcap$mle)
abline(v=coef(b_sig)["sigma1"],lty="dotted")

######## pf evaluation
p1_sig_pf <- apply(matrix(apply(ibpf_pf_sig_eval,2,mean),nrow=reps),2,max)
pf_prof_sig_range <- c(max(p1_sig_pf)-10,max(p1_sig_pf))
plot(y=p1_sig_pf,x=prof_sig_vals,ylim=pf_prof_sig_range,
  ylab="Log-likelihood",xlab="sigma1")
p1_sig_pf_mcap <- mcap(logLik=p1_sig_pf,parameter=prof_sig_vals)
lines(p1_sig_pf_mcap$fit$parameter,p1_sig_pf_mcap$fit$smoothed)

abline(v=p1_sig_pf_mcap$ci)
abline(v=p1_sig_pf_mcap$mle)
abline(v=coef(b_sig)["sigma1"],lty="dotted")





## ----model-he10-plot,eval=T,echo=F,fig.height=7, fig.width=6, out.width="6.5in", fig.cap = "Weekly measles case reports. (A) Data for four UK towns. (B) Simulated data."----
m_dir <- paste0("m_",i,"/")
if(!dir.exists(m_dir)) dir.create(m_dir)
stew(file=paste0(m_dir,"m-sim.rda"),{
  cat(capture.output(sessionInfo()),
    file=paste0(m_dir,"sessionInfo.txt"),sep="\n")
he10_model <- he10(U=4,dt=1/365,Tmax=switch(i,1955,1964),
  expandedParNames=c("R0"),
  towns_selected=c(1,2,11,12),
  basic_params = c(
    alpha =0.99,      iota=0,          R0=30,
    cohort=0.5,  amplitude=0.3,     gamma=52,
    sigma=52,           mu=0.02,  sigmaSE=0.05,
    rho=0.5,           psi=0.1,         g=800,
    S_0=0.036,         E_0=0.00007,   I_0=0.00006
  )
)
m <- simulate(he10_model,seed=27)
})
library(grid)
library(cowplot)
plot_data <- plot(he10_model,log=T)+xlab("")
plot_sim <- plot(m,log=T)+xlab("Year")
plot_grid(
  plot_data,
  plot_sim,
  labels=c("A","B"),
  nrow=2,ncol=1
)


## ----m-logLik-eval,cache=F,echo=F,eval=T--------------------------------------
stew(file=paste0(m_dir,"m-logLik.rda"),seed=10,{
  tic
  m_pf_logLik <- foreach(r=1:10,.combine=c) %dopar% {
    logLik(pfilter(m,switch(i,10,10000)))
  }
  m_bpf_logLik1 <- foreach(r=1:10,.combine=c) %dopar% {
    logLik(bpfilter(m,switch(i,10,10000),block_size=1))
  }
  m_bpf_logLik2 <- foreach(r=1:10,.combine=c) %dopar% {
    logLik(bpfilter(m,switch(i,10,10000),block_size=2))
  }
  m_bpf_logLik4 <- foreach(r=1:10,.combine=c) %dopar% {
    logLik(bpfilter(m,switch(i,10,10000),block_size=4))
  }
  toc
  m_lik_eval_time <- toc-tic
})



## ----m-prof-R0,eval=T,echo=F--------------------------------------------------
stew(file=paste0(m_dir,"m_prof_R02.rda"),seed=722,{
m_prof_R0_vals <- seq(from=25,to=35,length.out=switch(i,5,20))
tic <- Sys.time()
m_R0 <- m
m_R0_prof_u <- 1
reps <- switch(i,3,10)
param_R0_start_matrix <- matrix(coef(m_R0),byrow=TRUE,
  ncol=length(coef(m_R0)),nrow=reps*length(m_prof_R0_vals),
  dimnames=list(NULL,names(coef(m_R0))))
param_R0_start_matrix[,paste0("R0",m_R0_prof_u)] <- rep(m_prof_R0_vals,each=reps)
m_prof_R0_searches <- foreach(s=1:nrow(param_R0_start_matrix))%dopar%{
  rw_sd_call <- rep(list(0.02),times=U-1)
  names(rw_sd_call) <- lapply((1:U)[-m_R0_prof_u],function(x)paste0("R0",x))
  ibpf(m_R0,
    params=param_R0_start_matrix[s,],
    Nbpf=switch(i,2,50),Np=switch(i,10,2000),
    rw.sd=do.call(rw_sd,rw_sd_call),
    unitParNames="R0",
    sharedParNames=NULL,
    block_size=1,
    cooling.fraction.50=0.5
  )
}
toc <- Sys.time()
prof_R0_cores <- cores
})
m_prof_R0_time <- toc-tic


## ----m-prof-R0-lik,eval=T,echo=F----------------------------------------------
stew(file=paste0(m_dir,"m_prof_R0",m_R0_prof_u,"_lik.rda"),seed=868,{
tic <- Sys.time()
m_bpf_R0_eval <- foreach(mfit=m_prof_R0_searches,.combine=cbind) %dopar% {
  replicate(switch(i,2,10),
    logLik(bpfilter(mfit,block_size=1,Np=switch(i,10,5000)))
  )
}
m_pf_R0_eval <- foreach(mfit=m_prof_R0_searches,.combine=cbind) %dopar% {
  replicate(switch(i,2,10),
    logLik(pfilter(mfit,Np=switch(i,10,5000)))
  )
}
toc <- Sys.time()
prof_R0_eval_cores <- cores
})
m_prof_R0_eval_time <- toc-tic


## ----m-prof-R0-plot, fig.height=8, fig.width=6, out.width="5.5in", fig.cap = "Top: profile for $R0_1$ using an IBPF search with likelihood computed using BPF, for $K=4$ blocks each having 1 unit. Bottom: The same IBPF search with likelihood computed using the particle filter. Vertical lines show the MLE and a 95\\% confidence interval, with a dotted line at the true parameter value.",eval=T,echo=F----

par(mfrow=c(2,1))
par(mai=c(0.9,0.9,0.15,0.2))

##### window to avoid problematic particle filter behavior at extremes
prof_include <- (m_prof_R0_vals > 10) & (m_prof_R0_vals < 50)

##### bpf evaluation
m_p2_R0_bpf <- apply(matrix(apply(m_bpf_R0_eval,2,mean),nrow=reps),2,max)[prof_include]
bpf_prof_R0_range <- c(max(m_p2_R0_bpf)-80,max(m_p2_R0_bpf))
plot(y=m_p2_R0_bpf,x=m_prof_R0_vals[prof_include], ylim=bpf_prof_R0_range,
  ylab="Log-likelihood",xlab="")
m_p2_R0_bpf_mcap <- mcap(logLik=m_p2_R0_bpf,
  parameter=m_prof_R0_vals[prof_include])
lines(m_p2_R0_bpf_mcap$fit$parameter,m_p2_R0_bpf_mcap$fit$smoothed)

abline(v=m_p2_R0_bpf_mcap$ci)
abline(v=m_p2_R0_bpf_mcap$mle)
abline(v=coef(m_R0)[paste0("R0",m_R0_prof_u)],lty="dotted")

######## pf evaluation
m_p2_R0_pf <- apply(matrix(apply(m_pf_R0_eval,2,mean),nrow=reps),2,max)[prof_include]
pf_prof_R0_range <- c(max(m_p2_R0_pf)-80,max(m_p2_R0_pf))
plot(y=m_p2_R0_pf,x=m_prof_R0_vals[prof_include], ylim=pf_prof_R0_range,
  ylab="Log-likelihood",xlab="R0")
m_p2_R0_pf_mcap <- mcap(logLik=m_p2_R0_pf,parameter=m_prof_R0_vals[prof_include])
lines(m_p2_R0_pf_mcap$fit$parameter,m_p2_R0_pf_mcap$fit$smoothed)

abline(v=m_p2_R0_pf_mcap$ci)
abline(v=m_p2_R0_pf_mcap$mle)
abline(v=coef(m_R0)[paste0("R0",m_R0_prof_u)],lty="dotted")



## ----model-he10-g,eval=T,echo=F-----------------------------------------------
stew(file=paste0(m_dir,"m-g-sim.rda"),{
he10_model <- he10(U=4,dt=1/365,Tmax=switch(i,1955,1964),
  expandedParNames=c("g"),
  towns_selected=c(1,2,11,12),
  basic_params = c(
    alpha =0.99,      iota=0,          R0=30,
    cohort=0.5,  amplitude=0.3,     gamma=52,
    sigma=52,           mu=0.02,  sigmaSE=0.05,
    rho=0.5,           psi=0.1,         g=800,
    S_0=0.036,         E_0=0.00007,   I_0=0.00006
  )
)
m_g <- simulate(he10_model,seed=27)
if(any(dim(obs(m_g))!=dim(obs(m_R0))) || any(obs(m_g)-obs(m_R0) != 0)) stop("m_g and m_R0 should have the same synthetic data")
})


## ----m-slice-g,eval=T,echo=F--------------------------------------------------
stew(file=paste0(m_dir,"m_slice_g.rda"),seed=722,{
m_g_prof_u <- 4
slice_reps <- switch(i,2,10)
slice_points <- switch(i,5,20)
slice_g_vals <- rep(seq(from=100,to=2500,length.out=slice_points),each=slice_reps)
params_g_slice<- matrix(coef(m_g),
  nrow=slice_reps*slice_points,
  ncol=length(coef(m_g)),
  byrow=T,
  dimnames=list(NULL,names(coef(m_g))))
params_g_slice[,paste0("g",m_g_prof_u)] <- slice_g_vals
tic <- Sys.time()
m_slice_g_evals <- foreach(s=1:nrow(params_g_slice),.combine=c)%dopar%{
  logLik(bpfilter(m_g,
    params=params_g_slice[s,],
    Np=switch(i,10,5000),
    block_size=1
  ))
}
toc <- Sys.time()
m_slice_g_time <- toc-tic
m_slice_cores <- cores
})


## ----m_slice_g_plot,eval=T,echo=F,fig.height=4, fig.width=6, out.width="5.5in", fig.cap = "Slice for $g_4$ through the true parameter vector, using BPF evaluation with $K=4$ blocks each having 1 unit. Vertical lines show the MLE and a 95\\% confidence interval, with a dotted line at the true parameter value."----

plot(y=m_slice_g_evals,x=slice_g_vals,xlab="g",ylab="Log-likelihood")
m_slice_g_mcap <- mcap(logLik=m_slice_g_evals,parameter=slice_g_vals)
lines(m_slice_g_mcap$fit$parameter,m_slice_g_mcap$fit$smoothed)

abline(v=m_slice_g_mcap$ci)
abline(v=m_slice_g_mcap$mle)
abline(v=coef(m_g)[paste0("g",m_g_prof_u)],lty="dotted")



## ----diagnostics-setup,eval=T,echo=F------------------------------------------
diagnostics_dir <- paste0("diagnostics_",i,"/")
if(!dir.exists(diagnostics_dir)) dir.create(diagnostics_dir)


## ----pfilter-test-london,eval=T,echo=F,fig.height=6, fig.width=6, out.width="4.5in", fig.cap = "Comparison of effective sample size and conditional log-likelihood variance for measles in London"----
london <- he10(U=1)
london_diag <- bake(file=paste0(diagnostics_dir,"london.rds"),seed=230604, {
  pf_list <- foreach(rep=1:switch(i,6,100)) %dopar% pfilter(london,Np=switch(i,10,5000))
  ess=sapply(pf_list,function(x)x@eff.sample.size)
  ll=sapply(pf_list,function(x)x@cond.logLik)
  list(mess=apply(ess,1,mean),cllv=apply(ll,1,var),mll=apply(ll,1,mean))
})
london_ess <- london_diag$mess; london_cllv <- london_diag$cllv; londom_mll <- london_diag$mll
xy_range <- range(c(1/london_ess,london_cllv))
plot(x=london_cllv,y=1/london_ess,log="xy",ylab="1/ESS",xlab="Conditional log-likelihood variance",
  xlim=xy_range, ylim = xy_range)
abline(a=0,b=1)




## ----anomaly,echo=T,eval=F----------------------------------------------------
## benchmark <- arma_benchmark(he10_model)
## anomaly <- bcll - benchmark$cond


## ----bpfilter-diagnostics,eval=T,echo=F,fig.height=6.5, fig.width=6.5, out.width="6.5in", fig.cap = "Block particle filter diagnostics for London (A,C) and Birmingham (B,D)."----
bpf_diag <- bake(file=paste0(diagnostics_dir,"bpfilter.rds"),seed=230604, {
bpf_list <- foreach(rep=1:switch(i,6,40)) %dopar% {
  bpfilter(he10_model,Np=switch(i,10,2000),block_size=1)
}
bpf_ll  <- sapply(bpf_list,function(x)x@block.cond.loglik)
dim(bpf_ll) <- c(dim(bpf_list[[1]]@block.cond.loglik),length(bpf_list))
bcll <- apply(bpf_ll,c(1,2),mean)
bcllv <- apply(bpf_ll,c(1,2),var)
  list(bcll=bcll,bcllv=bcllv)
})
bcll <- bpf_diag$bcll; bcllv <- bpf_diag$bcllv
benchmark <- arma_benchmark(he10_model)
anomaly <- bcll - benchmark$cond
par(mfrow=c(2,2))
par(mai=c(1,0.8,0.1,0.1))
plot(x=bcllv[1,],y=anomaly[1,],log="x",ylab="anomaly",xlab="Conditional log-likelihood variance")
mtext("A",side=3,line=-1,adj=-0.24,cex=1.5)
plot(x=bcllv[2,],y=anomaly[2,],log="x",ylab="anomaly",xlab="Conditional log-likelihood variance")
mtext("B",side=3,line=-1,adj=-0.24,cex=1.5)
plot(y=anomaly[1,],x=time(he10_model),ylab="anomaly",xlab="Date")
mtext("C",side=3,line=-1,adj=-0.24,cex=1.5)
plot(y=anomaly[2,],x=time(he10_model),ylab="anomaly",xlab="Date")
mtext("D",side=3,line=-1,adj=-0.24,cex=1.5)


## ----traces,fig.height=6.5, fig.width=6.5, out.width="6.5in", fig.cap = "Trace plots for the correlated Gaussian random walk example."----
ibpf_traces <- pomp::melt(lapply(ibpf_mle_searches,pomp:::traces_internal))
ibpf_traces$iteration <- as.numeric(ibpf_traces$iteration)
ggplot(ibpf_traces, aes(x=iteration,y=value,group=.L1,color=factor(.L1))) +
  geom_line() +
  guides(color="none") +
  facet_wrap(~variable,scales="free_y")

